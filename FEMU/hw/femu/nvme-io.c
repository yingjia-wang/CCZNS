#include "./nvme.h"
#include "zns/zns.h"
#include "zns/znsftl.h"
#include "zstd.h"

#define CHECK_ZSTD(fn) \
    do { \
        size_t const err = (fn); \
        if (ZSTD_isError(err)) { \
            printf("CHECK failed: %s\n", ZSTD_getErrorName(err)); \
        } \
    } while (0) \

static uint16_t nvme_io_cmd(FemuCtrl *n, NvmeCmd *cmd, NvmeRequest *req);

static void nvme_update_sq_eventidx(const NvmeSQueue *sq)
{
    if (sq->eventidx_addr_hva) {
        *((uint32_t *)(sq->eventidx_addr_hva)) = sq->tail;
        return;
    }

    if (sq->eventidx_addr) {
        nvme_addr_write(sq->ctrl, sq->eventidx_addr, (void *)&sq->tail,
                        sizeof(sq->tail));
    }
}

static inline void nvme_copy_cmd(NvmeCmd *dst, NvmeCmd *src)
{
#if defined(__AVX__)
    __m256i *d256 = (__m256i *)dst;
    const __m256i *s256 = (const __m256i *)src;

    _mm256_store_si256(&d256[0], _mm256_load_si256(&s256[0]));
    _mm256_store_si256(&d256[1], _mm256_load_si256(&s256[1]));
#elif defined(__SSE2__)
    __m128i *d128 = (__m128i *)dst;
    const __m128i *s128 = (const __m128i *)src;

    _mm_store_si128(&d128[0], _mm_load_si128(&s128[0]));
    _mm_store_si128(&d128[1], _mm_load_si128(&s128[1]));
    _mm_store_si128(&d128[2], _mm_load_si128(&s128[2]));
    _mm_store_si128(&d128[3], _mm_load_si128(&s128[3]));
#else
    *dst = *src;
#endif
}

static void nvme_process_sq_io(void *opaque, int index_poller)
{
    NvmeSQueue *sq = opaque;
    FemuCtrl *n = sq->ctrl;

    hwaddr addr;
    NvmeCmd cmd;
    NvmeRequest *req;
    int processed = 0, rc;
    uint32_t next;

    nvme_update_sq_tail(sq);
    while (!(nvme_sq_empty(sq))) {
        if (sq->phys_contig) {
            addr = sq->dma_addr + sq->head * n->sqe_size;
            nvme_copy_cmd(&cmd, (void *)&(((NvmeCmd *)sq->dma_addr_hva)[sq->head]));
        } else {
            addr = nvme_discontig(sq->prp_list, sq->head, n->page_size,
                                  n->sqe_size);
            nvme_addr_read(n, addr, (void *)&cmd, sizeof(cmd));
        }
        nvme_inc_sq_head(sq);

        req = QTAILQ_FIRST(&sq->req_list);
        QTAILQ_REMOVE(&sq->req_list, req, entry);
        memset(&req->cqe, 0, sizeof(req->cqe));

        /* Coperd: record req->stime at earliest convenience */
        req->expire_time = req->stime = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        req->cqe.cid = cmd.cid;
        req->cmd_opcode = cmd.opcode;
        memcpy(&req->cmd, &cmd, sizeof(NvmeCmd));

        // poller transfers I/O to master, rather than actually process I/O before
        next = ++n->next_master;
        if (next < 1) next = 1;
        if (next > n->num_master) next = n->num_master;
        if (n->next_master >= n->num_master)
            n->next_master = 1;

        req->poller_id = index_poller;
        rc = femu_ring_enqueue(n->to_master[next], (void *)&req, 1);
        if (rc != 1)
            femu_err("enqueue failed, ret=%d\n", rc);

        processed++;
    }

    nvme_update_sq_eventidx(sq);
    sq->completed += processed;
}

static void nvme_post_cqe(NvmeCQueue *cq, NvmeRequest *req)
{
    FemuCtrl *n = cq->ctrl;
    NvmeSQueue *sq = req->sq;
    NvmeCqe *cqe = &req->cqe;
    uint8_t phase = cq->phase;
    hwaddr addr;

    if (n->print_log) {
        femu_debug("%s,req,lba:%lu,lat:%lu\n", n->devname, req->slba, req->reqlat);
    }
    cqe->status = cpu_to_le16((req->status << 1) | phase);
    cqe->sq_id = cpu_to_le16(sq->sqid);
    cqe->sq_head = cpu_to_le16(sq->head);

    if (cq->phys_contig) {
        addr = cq->dma_addr + cq->tail * n->cqe_size;
        ((NvmeCqe *)cq->dma_addr_hva)[cq->tail] = *cqe;
    } else {
        addr = nvme_discontig(cq->prp_list, cq->tail, n->page_size, n->cqe_size);
        nvme_addr_write(n, addr, (void *)cqe, sizeof(*cqe));
    }

    nvme_inc_cq_tail(cq);
}

static void nvme_process_cq_cpl(void *arg, int index_poller)
{
    FemuCtrl *n = (FemuCtrl *)arg;
    NvmeCQueue *cq = NULL;
    NvmeRequest *req = NULL;
    struct rte_ring *rp = n->to_ftl[index_poller];
    pqueue_t *pq = n->pq[index_poller];
    uint64_t now;
    int processed = 0;
    int rc;

    if (BBSSD(n)) {
        rp = n->to_poller[index_poller];
    }

    while (femu_ring_count(rp)) {
        req = NULL;
        rc = femu_ring_dequeue(rp, (void *)&req, 1);
        if (rc != 1) {
            femu_err("dequeue from to_poller request failed\n");
        }
        assert(req);

        pqueue_insert(pq, req);
    }

    while ((req = pqueue_peek(pq))) {
        now = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        if (now < req->expire_time) {
            break;
        }

        cq = n->cq[req->sq->sqid];
        if (!cq->is_active)
            continue;
        nvme_post_cqe(cq, req);
        QTAILQ_INSERT_TAIL(&req->sq->req_list, req, entry);
        pqueue_pop(pq);
        processed++;
        n->nr_tt_ios++;
        if (req->cmd_opcode == NVME_CMD_WRITE) n->nr_tt_write_ios++;
        else if (req->cmd_opcode == NVME_CMD_READ) n->nr_tt_read_ios++;
        if (now - req->expire_time >= 20000) {
            if (req->cmd_opcode == NVME_CMD_WRITE) n->nr_tt_write_late_ios++;
            else if (req->cmd_opcode == NVME_CMD_READ) n->nr_tt_read_late_ios++;
            femu_debug_time("actually %lu (expected %lu) cmd %d w %ld/%ld r %ld/%ld\n", 
                now - req->stime + req->process_time, req->reqlat + req->process_time, req->cmd_opcode,
                n->nr_tt_write_late_ios, n->nr_tt_write_ios,
                n->nr_tt_read_late_ios, n->nr_tt_read_ios);
        }
        n->should_isr[req->sq->sqid] = true;
    }

    if (processed == 0)
        return;

    switch (n->multipoller_enabled) {
    case 1:
        nvme_isr_notify_io(n->cq[index_poller]);
        break;
    default:
        for (int i = 1; i <= n->num_io_queues; i++) {
            if (n->should_isr[i]) {
                nvme_isr_notify_io(n->cq[i]);
                n->should_isr[i] = false;
            }
        }
        break;
    }
}

static void *nvme_poller(void *arg)
{
    FemuCtrl *n = ((NvmePollerThreadArgument *)arg)->n;
    int index = ((NvmePollerThreadArgument *)arg)->index;

    switch (n->multipoller_enabled) {
    case 1:
        while (1) {
            if ((!n->dataplane_started)) {
                usleep(1000);
                continue;
            }
            NvmeSQueue *sq = n->sq[index];
            NvmeCQueue *cq = n->cq[index];
            if (sq && sq->is_active && cq && cq->is_active) {
                nvme_process_sq_io(sq, index);
            }
            nvme_process_cq_cpl(n, index);
        }
        break;
    case 0:
        while (1) {
            if ((!n->dataplane_started)) {
                usleep(1000);
                continue;
            }

            for (int i = 1; i <= n->num_io_queues; i++) {
                NvmeSQueue *sq = n->sq[i];
                NvmeCQueue *cq = n->cq[i];
                if (sq && sq->is_active && cq && cq->is_active) {
                    nvme_process_sq_io(sq, index);
                }
            }
            nvme_process_cq_cpl(n, index);
        }
        break;
    default:
        femu_err("illegal multipoller_enabled param: %d\n", n->multipoller_enabled);
    }

    return NULL;
}

static uint16_t czns_map_dptr(FemuCtrl *n, size_t len, NvmeRequest *req)
{
    uint64_t prp1, prp2;

    switch (req->cmd.psdt) {
    case NVME_PSDT_PRP:
        prp1 = le64_to_cpu(req->cmd.dptr.prp1);
        prp2 = le64_to_cpu(req->cmd.dptr.prp2);

        return nvme_map_prp(&req->qsg, &req->iov, prp1, prp2, len, n);
    default:
        return NVME_INVALID_FIELD;
    }
}

static void czns_dma_load_data(FemuCtrl *n, size_t len, NvmeRequest *req, void *data)
{
    uint16_t status = czns_map_dptr(n, len, req);
    if (status)
        femu_err("czns_map_dptr failed\n");
    
    QEMUSGList *qsg = &req->qsg;
    int sg_cur_index = 0;
    dma_addr_t sg_cur_byte = 0;
    dma_addr_t cur_addr, cur_len;

    uint64_t data_oft = 0;

    while (sg_cur_index < qsg->nsg) {
        cur_addr = qsg->sg[sg_cur_index].base + sg_cur_byte;
        cur_len = qsg->sg[sg_cur_index].len - sg_cur_byte;

        if (dma_memory_rw(qsg->as, cur_addr, data + data_oft, cur_len, DMA_DIRECTION_TO_DEVICE, MEMTXATTRS_UNSPECIFIED)) {
            femu_err("czns dma error");
        }

        sg_cur_byte += cur_len;
        if (sg_cur_byte == qsg->sg[sg_cur_index].len) {
            sg_cur_byte = 0;
            ++sg_cur_index;
        }

        data_oft += cur_len;
    }

    qemu_sglist_destroy(qsg);
}

static uint16_t czns_map_mptr(FemuCtrl *n, size_t len, NvmeRequest *req)
{
    assert(req->cmd.psdt == NVME_PSDT_PRP);
    hwaddr mptr = le64_to_cpu(req->cmd.mptr);

    // We assume there is no cmb here!
    // assert(n->cmbsz == 0);

    pci_dma_sglist_init(&req->qsg, &n->parent_obj, 0);
    qemu_sglist_add(&req->qsg, mptr, len);

    return NVME_SUCCESS;
}

static void* czns_dma_load_metadata(FemuCtrl *n, size_t len, NvmeRequest *req)
{
    uint16_t status = czns_map_mptr(n, len, req);
    if (status)
        femu_err("czns_map_mptr failed\n");

    void *metadata = malloc(len); 
    backend_rw_data(metadata, &req->qsg, true);

    return metadata;
}

static void *nvme_master(void *arg)
{
    FemuCtrl *n = ((NvmePollerThreadArgument *)arg)->n;
    int index = ((NvmePollerThreadArgument *)arg)->index;

    NvmeNamespace *ns = NULL;
    NvmeRequest *req = NULL;
    NvmeRwCmd *rw = NULL;
    struct zonessd *ssd = n->znsssd;
    Zone *z = NULL;

    int rc, lb_size = ssd->sp.secsz, entry_cnt, first; // TO FIX
    uint32_t min_compressed_time = n->min_compressed_time, min_decompressed_time = n->min_decompressed_time;
    uint32_t compressed_time, decompressed_time;
    uint64_t slba, data_offset, pz_slba, max_data_length;
    uint32_t data_size, num_partitions, zone_idx, nlb, compressed_nlb, byte_offset, byte_increment, delta, *last_post = NULL;
    uint16_t *tmp = NULL, *tmp2 = NULL, intra_block_ofst; //, pad_sz;
    void *data = NULL, *metadata = NULL, *mb = NULL;
    chunk_info *ci = NULL;

    while (1) {
        if ((!n->dataplane_started)) {
            usleep(1000);
            continue;
        }

        if (!n->to_master[index] || !femu_ring_count(n->to_master[index]))
            continue;

        rc = femu_ring_dequeue(n->to_master[index], (void *)&req, 1);
        if (rc != 1) {
            femu_err("dequeue in master failed\n");
        }

        // identify compress/decompress flag in the command
        req->is_compress = req->is_decompress = (req->cmd.cdw12 & 0x10000);
        req->data_pages = 0;

        if (ZNSSD(n)) {      
            if (req->cmd.opcode == NVME_CMD_WRITE && req->is_compress) {
                ns = &n->namespaces[le32_to_cpu(req->cmd.nsid) - 1];
                rw = (NvmeRwCmd *)&req->cmd;
                nlb = (uint32_t)le16_to_cpu(rw->nlb) + 1;
                data_size = zns_l2b(ns, nlb);
                slba = le64_to_cpu(rw->slba);
                zone_idx = slba / n->zone_size;
                z = &(ssd->zone_array[zone_idx]);

                // update request start time
                compressed_time = (nlb <= 16? min_compressed_time : (min_compressed_time * nlb / 16));
                req->stime += compressed_time;
                req->process_time = compressed_time;

                // load data from DMA controller
                mb = n->mbe->logical_space;
                data_offset = zns_l2b(ns, n->zone_size * (n->num_zones + zone_idx * (n->space_amp - 1))) + z->byte_wp;

                femu_debug_time("begin dma transfer\n");
                data = mb + data_offset;
                czns_dma_load_data(n, data_size, req, data);

                metadata = czns_dma_load_metadata(n, n->oob_size * nlb, req);
                if (!metadata) 
                    femu_err("czns dma load metadata failed\n");   

                femu_debug_time("end dma transfer\n");

                tmp = (uint16_t *)metadata;
                num_partitions = *tmp;
                femu_debug_time("zone %u load %u bytes, has %u partitions\n", zone_idx, data_size, *tmp);
                if (!n->fixed_compression_ratio) {
                  tmp2 = malloc(n->oob_size * nlb);
                  memcpy(tmp2, tmp, n->oob_size * nlb); // save a copy of tmp

                  // let workers fetch compress tasks
                  pthread_spin_lock(&n->mst_mutex[index]);
                  n->mst_data[index] = data;
                  n->mst_metadata[index] = metadata;
                  n->mst_remaining[index] = num_partitions;
                  n->mst_finished[index] = 0;
                  n->mst_total_bytes[index] = data_size;
                  n->mst_start_byte[index] = 0;
                  pthread_spin_unlock(&n->mst_mutex[index]);

                  // harvest compressed data lengths from workers
                  while (true) {
                      if (n->mst_finished[index] == num_partitions)
                          break;
                      usleep(10);
                  }
                  femu_debug_time("ready to harvest\n");
                }

                // TO FIX: in logical blocks
                pz_slba = slba - zone_idx * ssd->sp.secs_per_pg * ssd->sp.pgs_per_zone;

                ci = z->cis[pz_slba];
                if (ci == NULL) {
                  femu_err("empty ci, check! (slba %lu zone_idx %u)\n", slba, zone_idx);
                }
                byte_offset = byte_increment = 0;

                if (!n->fixed_compression_ratio) {
                  tmp = (uint16_t *)metadata;
                  entry_cnt = 0;
                  for (int j = 1; j <= num_partitions; j++) {
                      // if (n->fixed_compression_ratio)
                      //   byte_increment = data_size / n->fixed_compression_ratio;
                      // else
                      byte_increment = tmp[j];

                      ci->lo[entry_cnt] = byte_offset % lb_size;
                      ci->po[entry_cnt] = z->byte_wp;
                      // if (n->fixed_compression_ratio)
                      //   ci->pre_sz[entry_cnt] = data_size;
                      // else
                      ci->pre_sz[entry_cnt] = tmp2[j];
                      ci->post_sz[entry_cnt] = byte_increment;
                      femu_debug_time("chunk[lo %u po %lu pre_sz %u post_sz %u]\n", 
                                    ci->lo[entry_cnt], ci->po[entry_cnt], 
                                    ci->pre_sz[entry_cnt], ci->post_sz[entry_cnt]);
                      ci->count ++;
                      last_post = &ci->post_sz[entry_cnt];
                      entry_cnt ++;

                      delta = (byte_offset + byte_increment) / lb_size - byte_offset / lb_size;
                      if (delta > 0) {
                        pz_slba += delta;
                        ci = z->cis[pz_slba];
                        if (ci == NULL)
                          femu_err("empty ci, check! (pz_slba %lu delta %u)\n", pz_slba, delta);
                        entry_cnt = 0;
                      }
                      if (entry_cnt % 5 == 0) { // TO FIX
                        // expand the space
                        ci->lo = realloc(ci->lo, sizeof(uint32_t) * (entry_cnt + 5));
                        ci->po = realloc(ci->po, sizeof(uint64_t) * (entry_cnt + 5));
                        ci->pre_sz = realloc(ci->pre_sz, sizeof(uint32_t) * (entry_cnt + 5));
                        ci->post_sz = realloc(ci->post_sz, sizeof(uint32_t) * (entry_cnt + 5));
                      }

                      // if (n->fixed_compression_ratio)
                      //   z->byte_wp += data_size;
                      // else
                      z->byte_wp += tmp2[j];

                      // update compressed length in CW command
                      tmp[j] = byte_increment;

                      byte_offset += byte_increment;
                      femu_debug_time("byte_increment %u byte_offset %u\n", byte_increment, byte_offset);
                  } 
                } else {
                  tmp = (uint16_t *)metadata;

                  byte_increment = data_size / n->fixed_compression_ratio;
                  ci->lo[0] = 0;
                  ci->po[0] = z->byte_wp;
                  ci->pre_sz[0] = data_size;
                  ci->post_sz[0] = byte_increment;

                  femu_debug_time("chunk[lo %u po %lu pre_sz %u post_sz %u]\n", 
                                ci->lo[0], ci->po[0], ci->pre_sz[0], ci->post_sz[0]);
                  ci->count ++;

                  pz_slba += byte_increment / lb_size;
                  ci = z->cis[pz_slba];
                  if (ci == NULL)
                    femu_err("empty ci, check! (pz_slba %lu delta %u)\n", pz_slba, delta);
                  z->byte_wp += data_size;
                  byte_offset += byte_increment;
                  femu_debug_time("byte_increment %u byte_offset %u\n", byte_increment, byte_offset);
                }

                if (n->page_padding && byte_offset % lb_size != 0) {
                  // pad data to integer logical blocks
                  delta = lb_size - (byte_offset % lb_size);
                  ssd->total_amp += delta;
                  tmp[num_partitions] += delta;
                  *last_post += delta;
                  byte_offset += delta;
                  femu_debug_time("lastly pad %u\n", delta);
                } else if (!n->page_padding) {
                  ;
                }

                req->metadata = metadata;
                req->data_pages = byte_offset / lb_size;
                tmp[num_partitions + 1] = byte_offset / lb_size;
                tmp[num_partitions + 2] = byte_offset % lb_size;
                free(tmp2);

                ssd->total_cwrite += req->data_pages * lb_size;

                femu_debug_time("compressed to %lu pages, padding? %d byte offset %u\n", 
                                req->data_pages, n->page_padding, byte_offset);         
            } 
            else if (req->cmd.opcode == NVME_CMD_READ && req->is_decompress) {
                // extract flags from the command
                bool flexible_mode = (req->cmd.cdw12 & 0x2000000);
                compressed_nlb = (req->cmd.cdw13 >> 8) & 0xff;
                intra_block_ofst = (req->cmd.cdw13 >> 16) & 0xffff;
                // pad_sz = req->cmd.cdw14 & 0xffff;

                ns = &n->namespaces[le32_to_cpu(req->cmd.nsid) - 1];
                rw = (NvmeRwCmd *)&req->cmd;
                nlb = (uint32_t)le16_to_cpu(rw->nlb) + 1;
                slba = le64_to_cpu(rw->slba);
                zone_idx = slba / n->zone_size;
                z = &(ssd->zone_array[zone_idx]);

                // update request start time
                decompressed_time = (compressed_nlb <= 16? min_decompressed_time : (min_decompressed_time * compressed_nlb / 16));
                req->stime += decompressed_time;
                req->process_time = decompressed_time;

                // load metadata from DMA controller
                femu_debug_time("begin dma transfer\n");
                metadata = czns_dma_load_metadata(n, n->oob_size * nlb, req);
                if (!metadata) 
                    femu_err("czns dma load metadata failed\n");      
                femu_debug_time("end dma transfer\n");

                tmp = (uint16_t *)metadata;
                num_partitions = 1;

                pz_slba = slba - zone_idx * ssd->sp.secs_per_pg * ssd->sp.pgs_per_zone;
                req->byte_length = 0;
                // if (flexible_mode)
                max_data_length = nlb * lb_size;
                // else
                //   max_data_length = compressed_nlb * lb_size - intra_block_ofst - pad_sz;
                femu_debug_time("dr zone_idx %u pz_slba %lu ibst %u flex? %d max %lu\n", zone_idx, pz_slba, 
                               intra_block_ofst, flexible_mode, max_data_length);

                first = true;
                for (uint64_t lba = pz_slba; lba <= pz_slba + compressed_nlb; lba++) {
                    ci = z->cis[lba];
                    if (ci == NULL)
                      femu_err("empty ci, check! (lba %lu)\n", lba);

                    for (int i = 0; i < ci->count; i++) {
                        femu_debug_time("loop ci from lba %lu lo %u pre %u post %u\n",
                                      lba, ci->lo[i], ci->pre_sz[i], ci->post_sz[i]);
                        if (flexible_mode && req->byte_length + ci->pre_sz[i] > max_data_length)
                          goto out;
                        if (!flexible_mode && !first)
                          goto out;
                        if (first && intra_block_ofst != ci->lo[i])
                          continue;
                        if (first) {
                          req->byte_start = ci->po[i];
                          first = false;
                          femu_debug_time("identify first chunk\n");
                        }
                        tmp[num_partitions++] = ci->pre_sz[i];
                        tmp[num_partitions++] = ci->post_sz[i];
                        if (ci->pre_sz[i] > 0xffff || ci->post_sz[i] > 0xffff)
                          femu_err("unexpectedly large chunk (pre: %u post: %u), check!\n", 
                                    ci->pre_sz[i], ci->post_sz[i]);
                        req->byte_length += ci->pre_sz[i];
                        femu_debug_time("identify chunk (pre: %u post: %u) from lba %lu i %d length %lu\n", 
                            ci->pre_sz[i], ci->post_sz[i], lba, i, req->byte_length);
                    }
                }
out:
                if (req->byte_length == 0) {
                  femu_err("no data fetched error, pz_slba %lu zone_idx %u cnlb %u\n", 
                                pz_slba, zone_idx, compressed_nlb);        
                }

                req->metadata = metadata;
                req->data_pages = compressed_nlb;
                *tmp = (num_partitions - 1) / 2;

                femu_debug_time("decompressed to %lu bytes\n", req->byte_length);
            }
        }

        // actually process I/O in master
        if (nvme_io_cmd(n, &(req->cmd), req) != NVME_SUCCESS)
            femu_err("I/O is not successful\n");
        femu_debug_time("after io processing\n");
        // send request back to cq poller
        rc = femu_ring_enqueue(n->to_ftl[req->poller_id], (void *)&req, 1);
        if (rc != 1) {
            femu_err("enqueue in master failed\n");
        }
    }

    return NULL;
}


static void *nvme_worker(void *arg)
{
    FemuCtrl *n = ((NvmePollerThreadArgument *)arg)->n;

    // TO FIX: each worker can accommodate at most 1MB data
    int i = 1, bsize = 1048576;
    void *tmp = malloc(bsize);
    uint16_t *md = NULL;
    uint32_t job_id, start_byte, num_bytes;
    size_t res;

    while (1) {
        if ((!n->dataplane_started)) {
            usleep(1000);
            continue;
        }

        // poll each master to fetch compress task
        if (n->mst_remaining[i] == 0) {
            if (++i > n->num_master) i = 1;
            continue;
        } else {
            pthread_spin_lock(&n->mst_mutex[i]);
            if (n->mst_remaining[i] == 0) {
                pthread_spin_unlock(&n->mst_mutex[i]);
                continue;
            }
            // get job info (exclusively between workers)
            md = (uint16_t *)n->mst_metadata[i];
            job_id = *md - n->mst_remaining[i] + 1;
            start_byte = n->mst_start_byte[i];
            num_bytes = n->fixed_compression_ratio? n->mst_total_bytes[i] : md[job_id]; // TO FIX
            n->mst_start_byte[i] = start_byte + num_bytes;
            n->mst_remaining[i] --;
            pthread_spin_unlock(&n->mst_mutex[i]);
        }

        // fixed_compression_ratio is for fio testing, where compressing large data by a single worker causes
        // the performance bottleneck.

        // int64_t time1 = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        if (n->fixed_compression_ratio) {
            res = num_bytes / n->fixed_compression_ratio;
        } else { // do compress data
            res = ZSTD_compress(tmp, bsize, n->mst_data[i] + start_byte, num_bytes, 1);
            CHECK_ZSTD(res);
        }
        // int64_t time2 = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        
        // write back the compressed length
        md[job_id] = res;
        pthread_spin_lock(&n->mst_mutex2[i]);
        n->mst_finished[i] ++;
        femu_debug_time("worker %d compress %u -> %d\n", ((NvmePollerThreadArgument *)arg)->index,
                      num_bytes, md[job_id]);
        pthread_spin_unlock(&n->mst_mutex2[i]);
    }

    return NULL;
}

static int cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)
{
    return (next > curr);
}

static pqueue_pri_t get_pri(void *a)
{
    return ((NvmeRequest *)a)->expire_time;
}

static void set_pri(void *a, pqueue_pri_t pri)
{
    ((NvmeRequest *)a)->expire_time = pri;
}

static size_t get_pos(void *a)
{
    return ((NvmeRequest *)a)->pos;
}

static void set_pos(void *a, size_t pos)
{
    ((NvmeRequest *)a)->pos = pos;
}

void nvme_create_poller(FemuCtrl *n)
{
    n->should_isr = g_malloc0(sizeof(bool) * (n->num_io_queues + 1));

    n->num_poller = n->multipoller_enabled ? n->num_io_queues : 1;
    /* Coperd: we put NvmeRequest into these rings */
    n->to_ftl = malloc(sizeof(struct rte_ring *) * (n->num_poller + 1));
    for (int i = 1; i <= n->num_poller; i++) {
        n->to_ftl[i] = femu_ring_create(FEMU_RING_TYPE_MP_SC, FEMU_MAX_INF_REQS);
        if (!n->to_ftl[i]) {
            femu_err("failed to create ring (n->to_ftl) ...\n");
            abort();
        }
        assert(rte_ring_empty(n->to_ftl[i]));
    }

    n->to_poller = malloc(sizeof(struct rte_ring *) * (n->num_poller + 1));
    for (int i = 1; i <= n->num_poller; i++) {
        n->to_poller[i] = femu_ring_create(FEMU_RING_TYPE_MP_SC, FEMU_MAX_INF_REQS);
        if (!n->to_poller[i]) {
            femu_err("failed to create ring (n->to_poller) ...\n");
            abort();
        }
        assert(rte_ring_empty(n->to_poller[i]));
    }

    n->pq = malloc(sizeof(pqueue_t *) * (n->num_poller + 1));
    for (int i = 1; i <= n->num_poller; i++) {
        n->pq[i] = pqueue_init(FEMU_MAX_INF_REQS, cmp_pri, get_pri, set_pri,
                               get_pos, set_pos);
        if (!n->pq[i]) {
            femu_err("failed to create pqueue (n->pq) ...\n");
            abort();
        }
    }

    n->to_master = malloc(sizeof(struct rte_ring *) * (n->num_master + 1));
    for (int i = 1; i <= n->num_master; i++) {
        n->to_master[i] = femu_ring_create(FEMU_RING_TYPE_MP_SC, FEMU_MAX_INF_REQS);
        if (!n->to_master[i]) {
            femu_err("failed to create ring (n->to_master) ...\n");
            abort();
        }
        assert(rte_ring_empty(n->to_master[i]));
    }

    n->poller = malloc(sizeof(QemuThread) * (n->num_poller + 1));
    NvmePollerThreadArgument *args = malloc(sizeof(NvmePollerThreadArgument) * (n->num_poller + 1));
    for (int i = 1; i <= n->num_poller; i++) {
        args[i].n = n;
        args[i].index = i;
        qemu_thread_create(&n->poller[i], "nvme-poller", nvme_poller, &args[i],
                           QEMU_THREAD_JOINABLE);
        femu_debug("nvme-poller [%d] created ...\n", i - 1);
    }
    free(args);

    n->master = malloc(sizeof(QemuThread) * (n->num_master + 1));
    args = malloc(sizeof(NvmePollerThreadArgument) * (n->num_master + 1));
    for (int i = 1; i <= n->num_master; i++) {
        args[i].n = n;
        args[i].index = i;
        qemu_thread_create(&n->master[i], "nvme-master", nvme_master, &args[i],
                           QEMU_THREAD_JOINABLE);
        femu_debug("nvme-master [%d] created ...\n", i - 1);
    }
    free(args);

    // Each master needs some data structures for dispatching and harvesting 
    // (de)compression jobs with high efficiency
    n->mst_data = malloc(sizeof(void *) * (n->num_master + 1));      
    n->mst_metadata = malloc(sizeof(void *) * (n->num_master + 1));
    n->mst_remaining = malloc(sizeof(uint32_t) * (n->num_master + 1));
    memset(n->mst_remaining, 0, sizeof(uint32_t) * (n->num_master + 1));
    n->mst_finished = malloc(sizeof(uint32_t) * (n->num_master + 1));
    n->mst_total_bytes = malloc(sizeof(uint32_t) * (n->num_master + 1));
    n->mst_start_byte = malloc(sizeof(uint32_t) * (n->num_master + 1));
    n->mst_mutex = malloc(sizeof(pthread_spinlock_t) * (n->num_master + 1));
    n->mst_mutex2 = malloc(sizeof(pthread_spinlock_t) * (n->num_master + 1));
    int ret;
    for (int i = 1; i <= n->num_master; i++) {
        n->mst_data[i] = NULL;
        n->mst_metadata[i] = NULL;
        ret = pthread_spin_init(&n->mst_mutex[i], PTHREAD_PROCESS_SHARED);
        if (ret != 0)
            femu_err("mst_mutex init error!\n");
        ret = pthread_spin_init(&n->mst_mutex2[i], PTHREAD_PROCESS_SHARED);
        if (ret != 0)
            femu_err("mst_mutex init error!\n");
    }

    n->worker = malloc(sizeof(QemuThread) * (n->num_worker + 1));
    args = malloc(sizeof(NvmePollerThreadArgument) * (n->num_worker + 1));
    for (int i = 1; i <= n->num_worker; i++) {
        args[i].n = n;
        args[i].index = i;
        qemu_thread_create(&n->worker[i], "nvme-worker", nvme_worker, &args[i],
                           QEMU_THREAD_JOINABLE);
        femu_debug("nvme-worker [%d] created ...\n", i - 1);
    }
    free(args);
}

uint16_t nvme_rw(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd, NvmeRequest *req)
{
    NvmeRwCmd *rw = (NvmeRwCmd *)cmd;
    uint16_t ctrl = le16_to_cpu(rw->control);
    uint32_t nlb  = le16_to_cpu(rw->nlb) + 1;
    uint64_t slba = le64_to_cpu(rw->slba);
    uint64_t prp1 = le64_to_cpu(rw->prp1);
    uint64_t prp2 = le64_to_cpu(rw->prp2);
    const uint8_t lba_index = NVME_ID_NS_FLBAS_INDEX(ns->id_ns.flbas);
    const uint16_t ms = le16_to_cpu(ns->id_ns.lbaf[lba_index].ms);
    const uint8_t data_shift = ns->id_ns.lbaf[lba_index].lbads;
    uint64_t data_size = (uint64_t)nlb << data_shift;
    uint64_t data_offset = slba << data_shift;
    uint64_t meta_size = nlb * ms;
    uint64_t elba = slba + nlb;
    uint16_t err;
    int ret;

    req->is_write = (rw->opcode == NVME_CMD_WRITE) ? 1 : 0;

    err = femu_nvme_rw_check_req(n, ns, cmd, req, slba, elba, nlb, ctrl,
                                 data_size, meta_size);
    if (err)
        return err;

    if (nvme_map_prp(&req->qsg, &req->iov, prp1, prp2, data_size, n)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                            offsetof(NvmeRwCmd, prp1), 0, ns->id);
        return NVME_INVALID_FIELD | NVME_DNR;
    }

    assert((nlb << data_shift) == req->qsg.size);

    req->slba = slba;
    req->status = NVME_SUCCESS;
    req->nlb = nlb;

    ret = backend_rw(n->mbe, &req->qsg, &data_offset, req->is_write);
    if (!ret) {
        return NVME_SUCCESS;
    }

    return NVME_DNR;
}

static uint16_t nvme_dsm(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                         NvmeRequest *req)
{
    uint32_t dw10 = le32_to_cpu(cmd->cdw10);
    uint32_t dw11 = le32_to_cpu(cmd->cdw11);
    uint64_t prp1 = le64_to_cpu(cmd->dptr.prp1);
    uint64_t prp2 = le64_to_cpu(cmd->dptr.prp2);

    if (dw11 & NVME_DSMGMT_AD) {
        uint16_t nr = (dw10 & 0xff) + 1;

        uint64_t slba;
        uint32_t nlb;
        NvmeDsmRange range[nr];

        if (dma_write_prp(n, (uint8_t *)range, sizeof(range), prp1, prp2)) {
            nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                                offsetof(NvmeCmd, dptr.prp1), 0, ns->id);
            return NVME_INVALID_FIELD | NVME_DNR;
        }

        req->status = NVME_SUCCESS;
        for (int i = 0; i < nr; i++) {
            slba = le64_to_cpu(range[i].slba);
            nlb = le32_to_cpu(range[i].nlb);
            if (slba + nlb > le64_to_cpu(ns->id_ns.nsze)) {
                nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_LBA_RANGE,
                                    offsetof(NvmeCmd, cdw10), slba + nlb, ns->id);
                return NVME_LBA_RANGE | NVME_DNR;
            }

            bitmap_clear(ns->util, slba, nlb);
        }
    }

    return NVME_SUCCESS;
}

static uint16_t nvme_compare(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                             NvmeRequest *req)
{
    NvmeRwCmd *rw = (NvmeRwCmd *)cmd;
    uint32_t nlb  = le16_to_cpu(rw->nlb) + 1;
    uint64_t slba = le64_to_cpu(rw->slba);
    uint64_t prp1 = le64_to_cpu(rw->prp1);
    uint64_t prp2 = le64_to_cpu(rw->prp2);

    uint64_t elba = slba + nlb;
    uint8_t lba_index = NVME_ID_NS_FLBAS_INDEX(ns->id_ns.flbas);
    uint8_t data_shift = ns->id_ns.lbaf[lba_index].lbads;
    uint64_t data_size = nlb << data_shift;
    uint64_t offset  = ns->start_block + (slba << data_shift);

    if ((slba + nlb) > le64_to_cpu(ns->id_ns.nsze)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_LBA_RANGE,
                            offsetof(NvmeRwCmd, nlb), elba, ns->id);
        return NVME_LBA_RANGE | NVME_DNR;
    }
    if (n->id_ctrl.mdts && data_size > n->page_size * (1 << n->id_ctrl.mdts)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                            offsetof(NvmeRwCmd, nlb), nlb, ns->id);
        return NVME_INVALID_FIELD | NVME_DNR;
    }
    if (nvme_map_prp(&req->qsg, &req->iov, prp1, prp2, data_size, n)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                            offsetof(NvmeRwCmd, prp1), 0, ns->id);
        return NVME_INVALID_FIELD | NVME_DNR;
    }
    if (find_next_bit(ns->uncorrectable, elba, slba) < elba) {
        return NVME_UNRECOVERED_READ;
    }

    for (int i = 0; i < req->qsg.nsg; i++) {
        uint32_t len = req->qsg.sg[i].len;
        uint8_t tmp[2][len];

        nvme_addr_read(n, req->qsg.sg[i].base, tmp[1], len);
        if (memcmp(tmp[0], tmp[1], len)) {
            qemu_sglist_destroy(&req->qsg);
            return NVME_CMP_FAILURE;
        }
        offset += len;
    }

    qemu_sglist_destroy(&req->qsg);

    return NVME_SUCCESS;
}

static uint16_t nvme_flush(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                           NvmeRequest *req)
{
    return NVME_SUCCESS;
}

static uint16_t nvme_write_zeros(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                                 NvmeRequest *req)
{
    NvmeRwCmd *rw = (NvmeRwCmd *)cmd;
    uint64_t slba = le64_to_cpu(rw->slba);
    uint32_t nlb  = le16_to_cpu(rw->nlb) + 1;

    if ((slba + nlb) > ns->id_ns.nsze) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_LBA_RANGE,
                            offsetof(NvmeRwCmd, nlb), slba + nlb, ns->id);
        return NVME_LBA_RANGE | NVME_DNR;
    }

    return NVME_SUCCESS;
}

static uint16_t nvme_write_uncor(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                                 NvmeRequest *req)
{
    NvmeRwCmd *rw = (NvmeRwCmd *)cmd;
    uint64_t slba = le64_to_cpu(rw->slba);
    uint32_t nlb  = le16_to_cpu(rw->nlb) + 1;

    if ((slba + nlb) > ns->id_ns.nsze) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_LBA_RANGE,
                            offsetof(NvmeRwCmd, nlb), slba + nlb, ns->id);
        return NVME_LBA_RANGE | NVME_DNR;
    }

    bitmap_set(ns->uncorrectable, slba, nlb);

    return NVME_SUCCESS;
}

static uint16_t nvme_io_cmd(FemuCtrl *n, NvmeCmd *cmd, NvmeRequest *req)
{
    NvmeNamespace *ns;
    uint32_t nsid = le32_to_cpu(cmd->nsid);

    if (nsid == 0 || nsid > n->num_namespaces) {
        femu_err("%s, NVME_INVALID_NSID %" PRIu32 "\n", __func__, nsid);
        return NVME_INVALID_NSID | NVME_DNR;
    }

    req->ns = ns = &n->namespaces[nsid - 1];

    switch (cmd->opcode) {
    case NVME_CMD_FLUSH:
        if (!n->id_ctrl.vwc || !n->features.volatile_wc) {
            return NVME_SUCCESS;
        }
        return nvme_flush(n, ns, cmd, req);
    case NVME_CMD_DSM:
        if (NVME_ONCS_DSM & n->oncs) {
            return nvme_dsm(n, ns, cmd, req);
        }
        return NVME_INVALID_OPCODE | NVME_DNR;
    case NVME_CMD_COMPARE:
        if (NVME_ONCS_COMPARE & n->oncs) {
            return nvme_compare(n, ns, cmd, req);
        }
        return NVME_INVALID_OPCODE | NVME_DNR;
    case NVME_CMD_WRITE_ZEROES:
        if (NVME_ONCS_WRITE_ZEROS & n->oncs) {
            return nvme_write_zeros(n, ns, cmd, req);
        }
        return NVME_INVALID_OPCODE | NVME_DNR;
    case NVME_CMD_WRITE_UNCOR:
        if (NVME_ONCS_WRITE_UNCORR & n->oncs) {
            return nvme_write_uncor(n, ns, cmd, req);
        }
        return NVME_INVALID_OPCODE | NVME_DNR;
    default:
        if (n->ext_ops.io_cmd) {
            return n->ext_ops.io_cmd(n, ns, cmd, req);
        }

        femu_err("%s, NVME_INVALID_OPCODE\n", __func__);
        return NVME_INVALID_OPCODE | NVME_DNR;
    }
}

void nvme_post_cqes_io(void *opaque)
{
    NvmeCQueue *cq = opaque;
    NvmeRequest *req, *next;
    int64_t cur_time, ntt = 0;
    int processed = 0;

    QTAILQ_FOREACH_SAFE(req, &cq->req_list, entry, next) {
        if (nvme_cq_full(cq)) {
            break;
        }

        cur_time = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
        if (cq->cqid != 0 && cur_time < req->expire_time) {
            ntt = req->expire_time;
            break;
        }

        nvme_post_cqe(cq, req);
        processed++;
    }

    if (ntt == 0) {
        ntt = qemu_clock_get_ns(QEMU_CLOCK_REALTIME) + CQ_POLLING_PERIOD_NS;
    }

    /* Only interrupt guest when we "do" complete some I/Os */
    if (processed > 0) {
        nvme_isr_notify_io(cq);
    }
}
