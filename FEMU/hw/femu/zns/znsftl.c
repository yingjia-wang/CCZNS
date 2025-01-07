#include "./znsftl.h"

static void zonessd_init_params(struct zonessd *ssd)
{
    struct zonessd_params *spp = &ssd->sp;

    
    spp->secsz = 4096;       // sector 4KB
    spp->secs_per_pg = 4;    // page 16KB
    spp->pgs_per_blk = 2048; // block 32MB
    spp->blks_per_pl = 64; 
    spp->pls_per_lun = 1;
    spp->luns_per_ch = 8;
    spp->nchs = 8;           // ssd 128GB

    spp->pg_rd_lat = NAND_READ_LATENCY;
    spp->pg_wr_lat = NAND_PROG_LATENCY;
    spp->blk_er_lat = NAND_ERASE_LATENCY;
    spp->ch_xfer_lat = 0;

    spp->pgs_per_pl = spp->pgs_per_blk * spp->blks_per_pl;
    spp->pgs_per_lun = spp->pgs_per_pl * spp->pls_per_lun;
    spp->pgs_per_ch = spp->pgs_per_lun * spp->luns_per_ch;
    spp->tt_pgs = spp->pgs_per_ch * spp->nchs;

    spp->blks_per_lun = spp->blks_per_pl * spp->pls_per_lun;
    spp->blks_per_ch = spp->blks_per_lun * spp->luns_per_ch;
    spp->tt_blks = spp->blks_per_ch * spp->nchs;

    spp->pls_per_ch =  spp->pls_per_lun * spp->luns_per_ch;
    spp->tt_pls = spp->pls_per_ch * spp->nchs;

    spp->tt_luns = spp->luns_per_ch * spp->nchs;

    /* line is special, put it at the end */
    spp->pgs_per_zone = ssd->zone_size_bs / spp->secsz / spp->secs_per_pg;
    spp->blks_per_zone = spp->pgs_per_zone / spp->pgs_per_blk;

    spp->enable_gc_delay = true;
}

static void zonessd_init_nand_blk(struct nand_block *blk, struct zonessd_params *spp)
{
    blk->npgs = spp->pgs_per_blk;
    blk->erase_cnt = 0;
    blk->wp = 0;
}

static void zonessd_init_nand_plane(struct nand_plane *pl, struct zonessd_params *spp)
{
    pl->nblks = spp->blks_per_pl;
    pl->blk = g_malloc0(sizeof(struct nand_block) * pl->nblks);
    for (int i = 0; i < pl->nblks; i++) {
        zonessd_init_nand_blk(&pl->blk[i], spp);
    }
}

static void zonessd_init_nand_lun(struct nand_lun *lun, struct zonessd_params *spp)
{
    lun->npls = spp->pls_per_lun;
    lun->pl = g_malloc0(sizeof(struct nand_plane) * lun->npls);
    for (int i = 0; i < lun->npls; i++) {
        zonessd_init_nand_plane(&lun->pl[i], spp);
    }
    lun->next_lun_avail_time = 0;
    lun->busy = false;
}

static void zonessd_init_ch(struct ssd_channel *ch, struct zonessd_params *spp)
{
    ch->nluns = spp->luns_per_ch;
    ch->lun = g_malloc0(sizeof(struct nand_lun) * ch->nluns);
    for (int i = 0; i < ch->nluns; i++) {
        zonessd_init_nand_lun(&ch->lun[i], spp);
    }
    ch->next_ch_avail_time = 0;
    ch->busy = false;
}

static void zonessd_init_zone_wp(Zone *Zone, int pos) 
{
    Zone->wp->zone = Zone;
    Zone->wp->pg = 0;
    Zone->wp->blk = 0;
}

static void zonessd_init_chunk_info(chunk_info *ci)
{
    ci->count = 0; 
    // TO FIX: reserve space for 5 entries by default, expand in need later
    int init_cnt = 5;
    ci->lo = malloc(sizeof(uint32_t) * init_cnt); // logical intra-block offset
    ci->po = malloc(sizeof(uint64_t) * init_cnt); // physical (expanded) offset 
    ci->pre_sz = malloc(sizeof(uint32_t) * init_cnt); // pre-compression size
    ci->post_sz = malloc(sizeof(uint32_t) * init_cnt); // post-compression size
}

static void zonessd_init_zone_array(struct zonessd *ssd) 
{
    struct zonessd_params *spp = &ssd->sp;

    ssd->zone_array = g_malloc0(sizeof(Zone) * ssd->zone_num);
    for (int i = 0; i < ssd->zone_num; i++) {
        Zone *z = &ssd->zone_array[i];
        z->wp = g_malloc0(sizeof(struct write_pointer));
        zonessd_init_zone_wp(&(ssd->zone_array[i]), i);
        z->pos = i;
        z->zone_size = ssd->zone_size;
        z->cis = g_malloc0(sizeof(chunk_info*) * spp->secs_per_pg * spp->pgs_per_zone);
        for (int j = 0; j < spp->secs_per_pg * spp->pgs_per_zone; j++) {
            z->cis[j] = g_malloc0(sizeof(struct chunk_info));
            zonessd_init_chunk_info(z->cis[j]);
        }
        z->buffer_persisted_time = g_malloc0(sizeof(uint64_t) * 10);
        z->buffer_cur_idx = 0;
        z->byte_wp = 0;
        z->lpn_wp = 0;
    }
}

static void zonessd_init_maptbl(struct zonessd *ssd)
{
    long long blk_cnt = 0;
    struct zonessd_params *spp = &ssd->sp;
    ssd->maptbl = g_malloc0(sizeof(struct zone_block) * ssd->zone_num);
    for (int zonenum = 0; zonenum < ssd->zone_num; zonenum++) {
        ssd->maptbl[zonenum].bl = g_malloc0(sizeof(struct block_locate) * ssd->sp.blks_per_zone);
        for (int blknum = 0; blknum < ssd->sp.blks_per_zone; blknum++) {
            // TO FIX
            // one-zone-all-dies configuration
            // ssd->maptbl[zonenum].bl[blknum].ch = (blk_cnt % spp->tt_pls) / spp->nchs;
            // ssd->maptbl[zonenum].bl[blknum].lun = (blk_cnt % spp->tt_pls % spp->nchs) / spp->pls_per_lun;
            
            // one-zone-half-dies configuration
            ssd->maptbl[zonenum].bl[blknum].ch = ((blk_cnt / spp->luns_per_ch) % (spp->nchs / 2));
            if (zonenum % 2 == 1) ssd->maptbl[zonenum].bl[blknum].ch += spp->nchs / 2;
            ssd->maptbl[zonenum].bl[blknum].lun = (blk_cnt % spp->luns_per_ch);
            femu_debug_time("zone %d blk %d ch %d lun %d\n", zonenum, blknum, 
                    ssd->maptbl[zonenum].bl[blknum].ch, 
                    ssd->maptbl[zonenum].bl[blknum].lun);
            blk_cnt++;
        }
    }
}

void zonessd_init(FemuCtrl *n) {
    struct zonessd *ssd = n->znsssd;
    struct zonessd_params *spp = &ssd->sp;

    assert(ssd);

    ssd->zone_size_bs = n->zone_size_bs;

    zonessd_init_params(ssd);

    /* initialize ssd internal layout architecture */
    ssd->ch = g_malloc0(sizeof(struct ssd_channel) * spp->nchs);
    for (int i = 0; i < spp->nchs; i++) {
        zonessd_init_ch(&ssd->ch[i], spp);
    }

    ssd->zone_size = n->zone_size / spp->secs_per_pg;
    ssd->zone_size_log2 = 0;
    if (is_power_of_2(n->zone_size)) {
        ssd->zone_size_log2 = 63 - clz64(ssd->zone_size);
    }

    ssd->zone_num = spp->tt_pgs / spp->pgs_per_zone;
    ssd->zone_array = g_malloc0(sizeof(Zone) * ssd->zone_num);

    ssd->total_amp = ssd->total_cwrite = ssd->total_write = 0;

    zonessd_init_zone_array(ssd);
    zonessd_init_maptbl(ssd);

    int ret = pthread_spin_init(&ssd->lock, PTHREAD_PROCESS_SHARED);
    if (ret != 0)
        femu_err("ssd->lock init error!\n");
}

static uint64_t ssd_advance_latency(struct zonessd *ssd, struct block_locate *bl, struct zone_cmd *ncmd) {
    int c = ncmd->cmd;
    uint64_t cmd_stime = (ncmd->stime == 0) ? \
        qemu_clock_get_ns(QEMU_CLOCK_REALTIME) : ncmd->stime;
    uint64_t nand_stime, lat = 0;

    struct zonessd_params *spp = &ssd->sp;
    struct nand_lun *lun = &ssd->ch[bl->ch].lun[bl->lun];

    femu_debug_time("hit ch %d lun %d\n", bl->ch, bl->lun);

    switch (c) {
    case NAND_READ:
        /* read: perform NAND cmd first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_rd_lat + DMA_TRANSFER_LATENCY + ECC_LATENCY;
        lat = lun->next_lun_avail_time - cmd_stime;
        break;

    case NAND_ERASE:
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->blk_er_lat;
        lat = lun->next_lun_avail_time - cmd_stime;
        break;

    default:
        femu_err("Unsupported NAND command: 0x%x\n", c);
    }
    
    return lat;
}

static uint64_t ssd_advance_stripe_latency(struct zonessd *ssd, struct block_locate *bl, struct zone_cmd *ncmd) {
    uint64_t cmd_stime = (ncmd->stime == 0) ? \
        qemu_clock_get_ns(QEMU_CLOCK_REALTIME) : ncmd->stime;
    uint64_t nand_stime, maxlat = 0, sublat;

    struct zonessd_params *spp = &ssd->sp;
    struct nand_lun *lun = NULL;

    femu_debug_time("hit ch [%d %d]\n", bl->ch - 3, bl->ch);

    // TO FIX
    for (int i = 3; i >= 0; i--) {
      for (int j = 0; j < 8; j++) {
        lun = &ssd->ch[bl->ch - i].lun[j];
        /* write: transfer data through channel first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                      lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat + DMA_TRANSFER_LATENCY;

        sublat = lun->next_lun_avail_time - cmd_stime;
        maxlat = (sublat > maxlat) ? sublat : maxlat;
      }
    }

    return maxlat;
}

static inline uint32_t zns_zone_idx(struct zonessd *ssd, uint64_t slba)
{
    return (ssd->zone_size_log2 > 0 ? slba >> ssd->zone_size_log2 : slba /
            ssd->zone_size);
}

static void check_addr(int num, int limit) {
    if (num >= limit) {
        printf("num : %d  limit : %d\n", num, limit);
    }
    assert(num >= 0 && num < limit);
}

static struct block_locate *zonessd_l2p(struct zonessd *ssd, uint64_t slba) 
{
    uint32_t zone_idx = zns_zone_idx(ssd, slba);
    uint32_t block_idx = (slba - zone_idx * ssd->sp.pgs_per_zone) % ssd->sp.blks_per_zone;
    check_addr(zone_idx, ssd->zone_num);
    struct block_locate *bl = &ssd->maptbl[zone_idx].bl[block_idx];

    return bl;
}

static void zonessd_clear_chunk_info(chunk_info *ci) {
    ci->count = 0;
    free(ci->lo);
    free(ci->po);
    free(ci->pre_sz);
    free(ci->post_sz);

    zonessd_init_chunk_info(ci);
}

void zonessd_reset(struct FemuCtrl *n, int pos, NvmeRequest *req) {
    struct zonessd *ssd = n->znsssd;
    Zone *z = &ssd->zone_array[pos];
    uint64_t maxlat = 0, sublat;

    pthread_spin_lock(&ssd->lock);
    for (int i = 0; i < ssd->sp.blks_per_zone; i++) {
        struct block_locate *bl = &ssd->maptbl[pos].bl[i];

        struct zone_cmd gce;
        gce.cmd = NAND_ERASE;
        gce.stime = req->stime;
        sublat = ssd_advance_latency(ssd, bl, &gce);
        maxlat = (sublat > maxlat) ? sublat : maxlat;
    }
    pthread_spin_unlock(&ssd->lock);

    z->wp->pg = 0;
    z->wp->blk = 0;

    for (int i = 0; i < ssd->sp.secs_per_pg * ssd->sp.pgs_per_zone; i++) {
      zonessd_clear_chunk_info(z->cis[i]);
    }
    z->byte_wp = 0;
    z->lpn_wp = 0;

    req->reqlat = maxlat;

    return ;
}

void zonessd_read(struct zonessd *ssd, NvmeRequest *req)
{
    struct zonessd_params *spp = &ssd->sp;
    NvmeRwCmd *rw = (NvmeRwCmd *)&req->cmd;
    uint64_t start_lpn = le64_to_cpu(rw->slba) / spp->secs_per_pg;
    uint32_t nlb = req->is_decompress? (req->data_pages - 1) : (uint32_t)le16_to_cpu(rw->nlb);
    uint64_t end_lpn = (le64_to_cpu(rw->slba) + nlb) / spp->secs_per_pg;

    femu_debug_time("read pages [%lu - %lu] slba %lu nlb %u\n", 
        start_lpn, end_lpn, le64_to_cpu(rw->slba), nlb);

    uint64_t maxlat = 0, sublat;

    pthread_spin_lock(&ssd->lock);
    for (uint64_t lpn = start_lpn; lpn <= end_lpn; lpn++) {
        struct block_locate *bl = zonessd_l2p(ssd, lpn);

        struct zone_cmd srd;
        srd.cmd = NAND_READ;
        srd.stime = req->stime;
        sublat = ssd_advance_latency(ssd, bl, &srd);
        maxlat = (sublat > maxlat) ? sublat : maxlat;
    }
    pthread_spin_unlock(&ssd->lock);

    req->reqlat = maxlat;
    req->expire_time = req->stime + maxlat;

    femu_debug_time("read lat %ld\n", req->reqlat);

    return ;
}

void zonessd_write(struct zonessd *ssd, NvmeRequest *req)
{
    struct zonessd_params *spp = &ssd->sp;
    NvmeRwCmd *rw = (NvmeRwCmd *)&req->cmd;
    uint64_t start_lpn = le64_to_cpu(rw->slba) / spp->secs_per_pg;
    uint32_t nlb = req->is_compress? (req->data_pages - 1) : (uint32_t)le16_to_cpu(rw->nlb);
    uint64_t end_lpn = (le64_to_cpu(rw->slba) + nlb) / spp->secs_per_pg;
    uint32_t zone_idx = zns_zone_idx(ssd, start_lpn);
    Zone *z = &ssd->zone_array[zone_idx];

    ssd->total_write += nlb * spp->secsz;

    femu_debug_time("write pages [%lu - %lu] zone %u slba %lu nlb %u\n", 
        start_lpn, end_lpn, zone_idx, le64_to_cpu(rw->slba), nlb);

    femu_debug_time("amp %lu/%lu/%lu\n", ssd->total_amp, ssd->total_cwrite, ssd->total_write);

    uint64_t maxlat = 0, sublat;

    pthread_spin_lock(&ssd->lock);
    for (uint64_t lpn = start_lpn; lpn <= end_lpn; lpn++) {
        struct block_locate *bl = zonessd_l2p(ssd, lpn);

        struct zone_cmd srd;
        srd.cmd = NAND_WRITE;
        srd.stime = req->stime;

        // indicate we can flush a stripe
        if (bl->lun == 7 && (bl->ch == 3 || bl->ch == 7) && z->lpn_wp < lpn) { // TO FIX
          sublat = ssd_advance_stripe_latency(ssd, bl, &srd);
          
          // record the buffer persisted time
          z->buffer_persisted_time[z->buffer_cur_idx] = req->stime + sublat;
          femu_debug_time("record bpt id %u time %lu\n", z->buffer_cur_idx, 
                        z->buffer_persisted_time[z->buffer_cur_idx]);
          if (++z->buffer_cur_idx == 10) z->buffer_cur_idx = 0; // TO FIX

          maxlat = (sublat > maxlat) ? sublat : maxlat;
        }
        z->lpn_wp = lpn;
    }
    pthread_spin_unlock(&ssd->lock);

    // TO FIX: two stripe buffers now
    int buffer_num = 1, idx = z->buffer_cur_idx;
    for (int i = 0; i < buffer_num; i++) 
      idx = idx == 0? 9 : (idx - 1);
    if (req->stime < z->buffer_persisted_time[idx]) {
      req->reqlat = z->buffer_persisted_time[idx] - req->stime;
      req->expire_time = req->stime + req->reqlat;
    } else {
      req->reqlat = 0;
      req->expire_time = req->stime;
    }

    femu_debug_time("check bpt id %u time %lu stime %lu write lat %lu\n", idx, 
                  z->buffer_persisted_time[idx], req->stime, req->reqlat);
    
    return ;
}

void zonessd_append(struct zonessd *ssd, NvmeRequest *req) {
    femu_err("FEMU does NOT support zone append now!\n");
    return ;
}
