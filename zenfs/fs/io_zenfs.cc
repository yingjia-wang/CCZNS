// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#if !defined(ROCKSDB_LITE) && !defined(OS_WIN)

#include "io_zenfs.h"

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <libzbd/zbd.h>
#include <linux/blkzoned.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "rocksdb/env.h"
#include "util/coding.h"

namespace ROCKSDB_NAMESPACE {

ZoneExtent::ZoneExtent(uint64_t start, uint64_t length, Zone* zone)
    : start_(start), length_(length), zone_(zone) {}

Status ZoneExtent::DecodeFrom(Slice* input) {
  if (input->size() != (sizeof(start_) + sizeof(length_)))
    return Status::Corruption("ZoneExtent", "Error: length missmatch");

  GetFixed64(input, &start_);
  GetFixed64(input, &length_);
  return Status::OK();
}

void ZoneExtent::EncodeTo(std::string* output) {
  PutFixed64(output, start_);
  PutFixed64(output, length_);
}

void ZoneExtent::EncodeJson(std::ostream& json_stream) {
  json_stream << "{";
  json_stream << "\"start\":" << start_ << ",";
  json_stream << "\"length\":" << length_;
  json_stream << "}";
}

enum ZoneFileTag : uint32_t {
  kFileID = 1,
  kFileNameDeprecated = 2,
  kFileSize = 3,
  kWriteLifeTimeHint = 4,
  kExtent = 5,
  kModificationTime = 6,
  kActiveExtentStart = 7,
  kIsSparse = 8,
  kLinkedFilename = 9,
};

void ZoneFile::EncodeTo(std::string* output, uint32_t extent_start) {
  PutFixed32(output, kFileID);
  PutFixed64(output, file_id_);

  PutFixed32(output, kFileSize);
  PutFixed64(output, file_size_);

  PutFixed32(output, kWriteLifeTimeHint);
  PutFixed32(output, (uint32_t)lifetime_);

  for (uint32_t i = extent_start; i < extents_.size(); i++) {
    std::string extent_str;

    PutFixed32(output, kExtent);
    extents_[i]->EncodeTo(&extent_str);
    PutLengthPrefixedSlice(output, Slice(extent_str));
  }

  PutFixed32(output, kModificationTime);
  PutFixed64(output, (uint64_t)m_time_);

  /* We store the current extent start - if there is a crash
   * we know that this file wrote the data starting from
   * active extent start up to the zone write pointer.
   * We don't need to store the active zone as we can look it up
   * from extent_start_ */
  PutFixed32(output, kActiveExtentStart);
  PutFixed64(output, extent_start_);

  if (is_sparse_) {
    PutFixed32(output, kIsSparse);
  }

  for (uint32_t i = 0; i < linkfiles_.size(); i++) {
    PutFixed32(output, kLinkedFilename);
    PutLengthPrefixedSlice(output, Slice(linkfiles_[i]));
  }
}

void ZoneFile::EncodeJson(std::ostream& json_stream) {
  json_stream << "{";
  json_stream << "\"id\":" << file_id_ << ",";
  json_stream << "\"size\":" << file_size_ << ",";
  json_stream << "\"hint\":" << lifetime_ << ",";
  json_stream << "\"extents\":[";

  for (const auto& name : GetLinkFiles())
    json_stream << "\"filename\":\"" << name << "\",";

  bool first_element = true;
  for (ZoneExtent* extent : extents_) {
    if (first_element) {
      first_element = false;
    } else {
      json_stream << ",";
    }
    extent->EncodeJson(json_stream);
  }
  json_stream << "]}";
}

Status ZoneFile::DecodeFrom(Slice* input) {
  uint32_t tag = 0;

  GetFixed32(input, &tag);
  if (tag != kFileID || !GetFixed64(input, &file_id_))
    return Status::Corruption("ZoneFile", "File ID missing");

  while (true) {
    Slice slice;
    ZoneExtent* extent;
    Status s;

    if (!GetFixed32(input, &tag)) break;

    switch (tag) {
      case kFileSize:
        if (!GetFixed64(input, &file_size_))
          return Status::Corruption("ZoneFile", "Missing file size");
        break;
      case kWriteLifeTimeHint:
        uint32_t lt;
        if (!GetFixed32(input, &lt))
          return Status::Corruption("ZoneFile", "Missing life time hint");
        lifetime_ = (Env::WriteLifeTimeHint)lt;
        break;
      case kExtent:
        extent = new ZoneExtent(0, 0, nullptr);
        GetLengthPrefixedSlice(input, &slice);
        s = extent->DecodeFrom(&slice);
        if (!s.ok()) {
          delete extent;
          return s;
        }
        extent->zone_ = zbd_->GetIOZone(extent->start_);
        if (!extent->zone_)
          return Status::Corruption("ZoneFile", "Invalid zone extent");
        extent->zone_->used_capacity_ += extent->length_;
        extents_.push_back(extent);
        break;
      case kModificationTime:
        uint64_t ct;
        if (!GetFixed64(input, &ct))
          return Status::Corruption("ZoneFile", "Missing creation time");
        m_time_ = (time_t)ct;
        break;
      case kActiveExtentStart:
        uint64_t es;
        if (!GetFixed64(input, &es))
          return Status::Corruption("ZoneFile", "Active extent start");
        extent_start_ = es;
        break;
      case kIsSparse:
        is_sparse_ = true;
        break;
      case kLinkedFilename:
        if (!GetLengthPrefixedSlice(input, &slice))
          return Status::Corruption("ZoneFile", "LinkFilename missing");

        if (slice.ToString().length() == 0)
          return Status::Corruption("ZoneFile", "Zero length Linkfilename");

        linkfiles_.push_back(slice.ToString());
        break;
      default:
        return Status::Corruption("ZoneFile", "Unexpected tag");
    }
  }

  MetadataSynced();
  return Status::OK();
}

Status ZoneFile::MergeUpdate(std::shared_ptr<ZoneFile> update, bool replace) {
  if (file_id_ != update->GetID())
    return Status::Corruption("ZoneFile update", "ID missmatch");

  SetFileSize(update->GetFileSize());
  SetWriteLifeTimeHint(update->GetWriteLifeTimeHint());
  SetFileModificationTime(update->GetFileModificationTime());

  if (replace) {
    ClearExtents();
  }

  std::vector<ZoneExtent*> update_extents = update->GetExtents();
  for (long unsigned int i = 0; i < update_extents.size(); i++) {
    ZoneExtent* extent = update_extents[i];
    Zone* zone = extent->zone_;
    zone->used_capacity_ += extent->length_;
    extents_.push_back(new ZoneExtent(extent->start_, extent->length_, zone));
  }
  extent_start_ = update->GetExtentStart();
  is_sparse_ = update->IsSparse();
  MetadataSynced();

  linkfiles_.clear();
  for (const auto& name : update->GetLinkFiles()) linkfiles_.push_back(name);

  return Status::OK();
}

ZoneFile::ZoneFile(ZonedBlockDevice* zbd, uint64_t file_id)
    : zbd_(zbd),
      active_zone_(NULL),
      extent_start_(NO_EXTENT),
      extent_filepos_(0),
      lifetime_(Env::WLTH_NOT_SET),
      io_type_(IOType::kUnknown),
      file_size_(0),
      file_id_(file_id),
      nr_synced_extents_(0),
      m_time_(0) {}

std::string ZoneFile::GetFilename() { return linkfiles_[0]; }
time_t ZoneFile::GetFileModificationTime() { return m_time_; }

uint64_t ZoneFile::GetFileSize() { return file_size_; }
void ZoneFile::SetFileSize(uint64_t sz) { file_size_ = sz; }
void ZoneFile::SetFileModificationTime(time_t mt) { m_time_ = mt; }
void ZoneFile::SetIOType(IOType io_type) { io_type_ = io_type; }

ZoneFile::~ZoneFile() {
  ClearExtents();
  IOStatus s = CloseWR();
  if (!s.ok()) {
    zbd_->SetZoneDeferredStatus(s);
  }
}

void ZoneFile::ClearExtents() {
  for (auto e = std::begin(extents_); e != std::end(extents_); ++e) {
    Zone* zone = (*e)->zone_;

    assert(zone && zone->used_capacity_ >= (*e)->length_);
    zone->used_capacity_ -= (*e)->length_;
    delete *e;
  }
  extents_.clear();
}

IOStatus ZoneFile::CloseWR() {
  IOStatus s;
  if (open_for_wr_) {
    /* Mark up the file as being closed */
    extent_start_ = NO_EXTENT;
    s = PersistMetadata();
    if (!s.ok()) return s;
    open_for_wr_ = false;
    s = CloseActiveZone();
  }
  return s;
}

IOStatus ZoneFile::CloseActiveZone() {
  IOStatus s = IOStatus::OK();
  if (active_zone_) {
    bool full = active_zone_->IsFull();
    s = active_zone_->Close();
    ReleaseActiveZone();
    if (!s.ok()) {
      return s;
    }
    zbd_->PutOpenIOZoneToken();
    if (full) {
      zbd_->PutActiveIOZoneToken();
    }
  }
  return s;
}

void ZoneFile::OpenWR(MetadataWriter* metadata_writer) {
  open_for_wr_ = true;
  metadata_writer_ = metadata_writer;
}

bool ZoneFile::IsOpenForWR() { return open_for_wr_; }

IOStatus ZoneFile::PersistMetadata() {
  assert(metadata_writer_ != NULL);
  return metadata_writer_->Persist(this);
}

ZoneExtent* ZoneFile::GetExtent(uint64_t file_offset, uint64_t* dev_offset) {
  for (unsigned int i = 0; i < extents_.size(); i++) {
    if (file_offset < extents_[i]->length_) {
      *dev_offset = extents_[i]->start_ + file_offset;
      return extents_[i];
    } else {
      file_offset -= extents_[i]->length_;
    }
  }
  return NULL;
}

IOStatus ZoneFile::PositionedRead(uint64_t offset, size_t n, Slice* result,
                                  char* scratch, const IOOptions& opts) {
  ZenFSMetricsLatencyGuard guard(zbd_->GetMetrics(), ZENFS_READ_LATENCY,
                                 Env::Default());
  zbd_->GetMetrics()->ReportQPS(ZENFS_READ_QPS, 1);

  ReadLock lck(this);

  uint64_t* alens_addr = opts.alens_addr;
  bool decompressed_read = opts.decompressed_read;
  bool flexible_mode = opts.flexible_mode;
  uint64_t actual_size = opts.actual_size;
#ifndef NDBG
  printf("in ZoneFile::PositionedRead, o %lu n %lu dr? %d fx? %d as %lu fs %lu fn %s\n", offset, n, 
          decompressed_read, flexible_mode, actual_size, GetFileSize(), GetFilename().c_str());
#endif
  // int f = zbd_->GetReadFD();
  int f_direct = zbd_->GetReadDirectFD();
  char* ptr;
  uint64_t r_off;
  size_t r_sz;
  ssize_t r = 0;
  size_t read = 0, aread = 0;
  ZoneExtent* extent;
  uint64_t extent_end;
  IOStatus s;

  if (offset >= file_size_) {
    if (offset > file_size_) printf("error! offset beyond file_size_\n");
    *result = Slice(scratch, 0);
    return IOStatus::OK();
  }

  r_off = 0;
  extent = GetExtent(offset, &r_off);
  if (!extent) {
    printf("error! read start beyond end of (synced) file data\n");
    *result = Slice(scratch, 0);
    return s;
  }
  extent_end = extent->start_ + extent->length_;

  /* Limit read size to end of file */
  if ((offset + n) > file_size_)
    r_sz = file_size_ - offset;
  else
    r_sz = n;

  ptr = scratch;

  void *metadata = nullptr;
  uint16_t *tmp = nullptr;
  size_t upper = 1048576, tmp_len = 0, read_time = 0; // TO FIX
  uint32_t block_sz = zbd_->GetBlockSize();

  while (read != r_sz) {
    read_time ++;
    size_t pread_sz = r_sz - read;

    if (pread_sz > upper) {
      printf("check! unexpected size %lu bytes!\n", pread_sz);
      // printf("in ZoneFile::PositionedRead, o %lu n %lu dr? %d fx? %d as %lu fs %lu fn %s\n", offset, n, 
      //        decompressed_read, flexible_mode, actual_size, GetFileSize(), GetFilename().c_str());     
      pread_sz = upper;
    }

    if ((pread_sz + r_off) > extent_end) {
#ifndef NDBG
      printf("(pread_sz + r_off) beyond extent_end, pread_sz %lu r_off %lu extent_end %lu\n",
              pread_sz, r_off, extent_end);
#endif
      pread_sz = extent_end - r_off;
    }

    /* Align read request to the blocksize boundaries */
    size_t pread_sz_aligned = pread_sz;
    pread_sz_aligned += r_off % block_sz;
    if ((r_off + pread_sz) % block_sz != 0)
      pread_sz_aligned += block_sz - (r_off + pread_sz) % block_sz;

    // preserve enough space since data length AFTER decompressed is often
    // larger than that BEFORE decompressed
    uint32_t data_len = pread_sz_aligned;
    uint32_t metadata_len = 0, oob_sz = zbd_->GetOOBSize();

    if (decompressed_read) {
      if (flexible_mode) {
        // Sometimes the estimation of data length even can not accommodate
        // one data block (if it is large), so we should additionally determine
        // based on the actual size (which is exactly of the first data block).
        data_len = std::max(pread_sz * 3, actual_size << 1);
        // TO FIX: configure a minimum size to prevent overflow
        data_len = std::max(data_len, (uint32_t)(0x4000));
      } else {
        data_len = actual_size;
      }
      if (data_len % block_sz != 0)
        data_len += block_sz - data_len % block_sz;
      metadata_len = data_len / block_sz * oob_sz;
      metadata = malloc(metadata_len);
    }

    uint64_t r_off_page = r_off >> 12;
    uint64_t intra_block_ofst = r_off % block_sz;
    uint32_t cdw12, cdw13, cdw14;
    uint32_t ori_nlb = (pread_sz_aligned >> 12);
    if (decompressed_read) {
      cdw12 = ((data_len >> 12) - 1) | (1 << 16);
      if (flexible_mode) cdw12 |= (1 << 25);
      cdw13 = (ori_nlb << 8) | (intra_block_ofst << 16);
      cdw14 = flexible_mode? 0 : (data_len - actual_size);
    } else {
      cdw12 = ori_nlb - 1;
      cdw13 = 0;
      cdw14 = 0;
    }

    struct nvme_passthru_cmd cmd = {
      .opcode = 2, // read
      .nsid = 1,
      .metadata = (__u64)metadata,
      .addr = (__u64)ptr,
      .metadata_len = metadata_len,
      .data_len = data_len,
      .cdw10 = (__u32)(r_off_page & 0xffffffff),
      .cdw11 = (__u32)(r_off_page >> 32),
      .cdw12 = cdw12,
      .cdw13 = cdw13,
      .cdw14 = cdw14,
    };

    int ret;
    // r = pread(f, ptr, pread_sz, r_off);
    // printf("read! r_off %lu r_off_page %lu prsz %lu prsz_a %lu ib_ofst %lu w10 %x w11 %x w12 %x w13 %x w14 %x\n", r_off, r_off_page, 
    //        pread_sz, pread_sz_aligned, intra_block_ofst, cmd.cdw10, cmd.cdw11, cmd.cdw12, cmd.cdw13, cmd.cdw14);

    // auto start = std::chrono::high_resolution_clock::now(); 
    // zbd_->ioctl_.lock();
    ret = ioctl(f_direct, NVME_IOCTL_IO_CMD, &cmd);
    // zbd_->ioctl_.unlock();
    // auto end = std::chrono::high_resolution_clock::now(); 
    // std::chrono::duration<double, std::milli> d = end - start;
    // zbd_->avg_read_io_time = (zbd_->avg_read_io_time * zbd_->read_io_count + d.count()) / (zbd_->read_io_count + 1);
    // zbd_->read_io_count ++;

    if (ret != 0) {
      if (ret == -1 && errno == EINTR)
        continue;
      r = ret;
      printf("read error(%d %d)! r_off %lu -> %lu prsz %lu len %u ori %u\n", ret, errno, 
              r_off, r_off_page, pread_sz, data_len, ori_nlb);
      break;
    }

    uint16_t *md = (uint16_t *)metadata;
    uint32_t asize = 0;
    if (decompressed_read) {
      for (int i = 1; i <= *md; i++)
        asize += md[(i << 1) - 1];
      if ((asize > data_len) || ((uint32_t)(*md) * 4 + 2 > metadata_len))
        printf("space overflow, asize %u data_len %u fx? %d num %d mlen %u\n", asize, data_len, 
                flexible_mode, *md, metadata_len);
    }

    ptr += decompressed_read? asize : pread_sz;
    aread += decompressed_read? asize : pread_sz;
    read += pread_sz;
    r_off += pread_sz;
    
    if (read != r_sz && r_off == extent_end) {
      extent = GetExtent(offset + read, &r_off);
      if (!extent) {
        printf("error! read beyond end of (synced) file data again\n");
        break;
      }
      r_off = extent->start_;
      extent_end = extent->start_ + extent->length_;
#ifndef NDBG
      printf("get a new extent: r_off %lu extent_end %lu\n", r_off, extent_end);
#endif
    }

    if (alens_addr) {
      if (read == r_sz && read_time == 1) {
        tmp = (uint16_t*)metadata;
      } else if (read != r_sz && read_time == 1) {
        // read I/O crosses different extents, thus we need to concatenate metadata together;
        // this condition happens during the first read I/O
        tmp = (uint16_t*)malloc(metadata_len);
        memcpy(tmp, md, ((*md << 1) + 1) * sizeof(uint16_t));
        tmp_len += ((*md << 1) + 1);
        free(metadata);
      } else {
        // read I/O crosses different extents, thus we need to concatenate metadata together;
        tmp = (uint16_t*)realloc(tmp, tmp_len + metadata_len);
        *tmp = *tmp + *md;
        memcpy(tmp + tmp_len, md + 1, (*md << 1) * sizeof(uint16_t));
        tmp_len += (*md << 1);
        free(metadata);
      }
    }
  }

  if (r != 0) {
    s = IOStatus::IOError("pread error\n");
    read = 0;
  }

  if (alens_addr) {
    *alens_addr = (__u64)tmp;
  }

#ifndef NDBG
  if(decompressed_read) 
    printf("got aread %lu vs read %lu\n", aread, read);
#endif
  // if (zbd_->read_io_count % 1000000000)
  //   printf("read epoch %lu: avg %lf ms\n", zbd_->read_io_count, zbd_->avg_read_io_time);
  
  *result = Slice((char*)scratch, aread);
  return s;
}

void ZoneFile::PushExtent() {
  uint64_t length;

  assert(file_size_ >= extent_filepos_);

  if (!active_zone_) return;

  length = file_size_ - extent_filepos_;
  if (length == 0) return;

  assert(length <= (active_zone_->wp_ - extent_start_));
  extents_.push_back(new ZoneExtent(extent_start_, length, active_zone_));

  active_zone_->used_capacity_ += length;
  extent_start_ = active_zone_->wp_;
  extent_filepos_ = file_size_;
}

IOStatus ZoneFile::AllocateNewZone() {
  Zone* zone;
  IOStatus s = zbd_->AllocateIOZone(lifetime_, io_type_, &zone);

  if (!s.ok()) return s;
  if (!zone) {
    return IOStatus::NoSpace("Zone allocation failure\n");
  }
  SetActiveZone(zone);
  extent_start_ = active_zone_->wp_;
  extent_filepos_ = file_size_;

  /* Persist metadata so we can recover the active extent using
     the zone write pointer in case there is a crash before syncing */
  return PersistMetadata();
}

/* Byte-aligned writes without a sparse header */
IOStatus ZoneFile::BufferedAppend(char* buffer, uint32_t data_size) {
  uint32_t left = data_size;
  uint32_t wr_size;
  uint32_t block_sz = GetBlockSize();
  IOStatus s;
  // printf("in ZoneFile::BufferedAppend, size %u fn %s\n", data_size, GetFilename().c_str());
  if (active_zone_ == NULL) {
    s = AllocateNewZone();
    if (!s.ok()) return s;
  }

  while (left) {
    wr_size = left;
    if (wr_size > active_zone_->capacity_) wr_size = active_zone_->capacity_;

    /* Pad to the next block boundary if needed */
    uint32_t align = wr_size % block_sz;
    uint32_t pad_sz = 0;

    if (align) pad_sz = block_sz - align;

    /* the buffer size s aligned on block size, so this is ok*/
    if (pad_sz) memset(buffer + wr_size, 0x0, pad_sz);

    uint64_t extent_length = wr_size;

    s = active_zone_->Append(buffer, wr_size + pad_sz);
    if (!s.ok()) return s;

    extents_.push_back(
        new ZoneExtent(extent_start_, extent_length, active_zone_));

    extent_start_ = active_zone_->wp_;
    active_zone_->used_capacity_ += extent_length;
    file_size_ += extent_length;
    left -= extent_length;

    if (active_zone_->capacity_ == 0) {
      s = CloseActiveZone();
      if (!s.ok()) {
        return s;
      }
      if (left) {
        memmove((void*)(buffer), (void*)(buffer + wr_size), left);
      }
      s = AllocateNewZone();
      if (!s.ok()) return s;
    }
  }

  return IOStatus::OK();
}

/* Byte-aligned, sparse writes with inline metadata
   the caller reserves 8 bytes of data for a size header */
IOStatus ZoneFile::SparseAppend(char* sparse_buffer, uint32_t data_size) {
  uint32_t left = data_size;
  uint32_t wr_size;
  uint32_t block_sz = GetBlockSize();
  IOStatus s;
  // printf("in ZoneFile::SparseAppend, size %u fn %s\n", data_size, GetFilename().c_str());
  if (active_zone_ == NULL) {
    s = AllocateNewZone();
    if (!s.ok()) return s;
  }

  while (left) {
    wr_size = left + ZoneFile::SPARSE_HEADER_SIZE;
    if (wr_size > active_zone_->capacity_) wr_size = active_zone_->capacity_;

    /* Pad to the next block boundary if needed */
    uint32_t align = wr_size % block_sz;
    uint32_t pad_sz = 0;

    if (align) pad_sz = block_sz - align;

    /* the sparse buffer has block_sz extra bytes tail allocated for padding, so
     * this is safe */
    if (pad_sz) memset(sparse_buffer + wr_size, 0x0, pad_sz);

    uint64_t extent_length = wr_size - ZoneFile::SPARSE_HEADER_SIZE;
    EncodeFixed64(sparse_buffer, extent_length);

    s = active_zone_->Append(sparse_buffer, wr_size + pad_sz);
    if (!s.ok()) return s;

    extents_.push_back(
        new ZoneExtent(extent_start_ + ZoneFile::SPARSE_HEADER_SIZE,
                       extent_length, active_zone_));

    extent_start_ = active_zone_->wp_;
    active_zone_->used_capacity_ += extent_length;
    file_size_ += extent_length;
    left -= extent_length;

    if (active_zone_->capacity_ == 0) {
      s = CloseActiveZone();
      if (!s.ok()) {
        return s;
      }
      if (left) {
        memmove((void*)(sparse_buffer + ZoneFile::SPARSE_HEADER_SIZE),
                (void*)(sparse_buffer + wr_size), left);
      }
      s = AllocateNewZone();
      if (!s.ok()) return s;
    }
  }

  return IOStatus::OK();
}

/* Assumes that data and size are block aligned */
IOStatus ZoneFile::Append(void* data, int data_size, void *metadata, int metadata_num, size_t* asize) {
  uint32_t left = data_size;
  uint32_t wr_size, offset = 0;
  IOStatus s = IOStatus::OK();

  if (!active_zone_) {
    s = AllocateNewZone();
    if (!s.ok()) return s;
  }

  while (left) {
    if (active_zone_->capacity_ == 0) {
      PushExtent();

      s = CloseActiveZone();
      if (!s.ok()) {
        return s;
      }

      s = AllocateNewZone();
      if (!s.ok()) return s;
    }

    // It is possible that the original data is NOT able to store in active_zone_,
    // but the compressed data can.
    // Let's simplify this situation, if the original data can NOT be stored in a
    // zone, this zone will be closed and a new active zone will be allocated.
    // Original codes:
    //     wr_size = left;
    //     if (wr_size > active_zone_->capacity_) wr_size = active_zone_->capacity_;
    if (left > active_zone_->capacity_) {
      wr_size = 0;
      // active_zone_->used_capacity_ += active_zone_->capacity_;
      active_zone_->capacity_ = 0;
    } else {
      wr_size = left;

      s = active_zone_->Append((char*)data + offset, left, metadata, metadata_num, asize);
      if (!s.ok()) return s;

      file_size_ += *asize;
      left -= wr_size;
      offset += wr_size;
    }
  }

  return IOStatus::OK();
}

IOStatus ZoneFile::RecoverSparseExtents(uint64_t start, uint64_t end,
                                        Zone* zone) {
  /* Sparse writes, we need to recover each individual segment */
  IOStatus s;
  uint32_t block_sz = GetBlockSize();
  int f = zbd_->GetReadFD();
  uint64_t next_extent_start = start;
  char* buffer;
  int recovered_segments = 0;
  int ret;

  ret = posix_memalign((void**)&buffer, sysconf(_SC_PAGESIZE), block_sz);
  if (ret) {
    return IOStatus::IOError("Out of memory while recovering");
  }

  while (next_extent_start < end) {
    uint64_t extent_length;

    ret = pread(f, (void*)buffer, block_sz, next_extent_start);
    if (ret != (int)block_sz) {
      s = IOStatus::IOError("Unexpected read error while recovering");
      break;
    }

    extent_length = DecodeFixed64(buffer);
    if (extent_length == 0) {
      s = IOStatus::IOError("Unexexpeted extent length while recovering");
      break;
    }
    recovered_segments++;

    zone->used_capacity_ += extent_length;
    extents_.push_back(new ZoneExtent(next_extent_start + SPARSE_HEADER_SIZE,
                                      extent_length, zone));

    uint64_t extent_blocks = (extent_length + SPARSE_HEADER_SIZE) / block_sz;
    if ((extent_length + SPARSE_HEADER_SIZE) % block_sz) {
      extent_blocks++;
    }
    next_extent_start += extent_blocks * block_sz;
  }

  free(buffer);
  return s;
}

IOStatus ZoneFile::Recover() {
  /* If there is no active extent, the file was either closed gracefully
     or there were no writes prior to a crash. All good.*/
  if (!HasActiveExtent()) return IOStatus::OK();

  /* Figure out which zone we were writing to */
  Zone* zone = zbd_->GetIOZone(extent_start_);

  if (zone == nullptr) {
    return IOStatus::IOError(
        "Could not find zone for extent start while recovering");
  }

  if (zone->wp_ < extent_start_) {
    return IOStatus::IOError("Zone wp is smaller than active extent start");
  }

  /* How much data do we need to recover? */
  uint64_t to_recover = zone->wp_ - extent_start_;

  /* Do we actually have any data to recover? */
  if (to_recover == 0) {
    /* Mark up the file as having no missing extents */
    extent_start_ = NO_EXTENT;
    return IOStatus::OK();
  }

  /* Is the data sparse or was it writted direct? */
  if (is_sparse_) {
    IOStatus s = RecoverSparseExtents(extent_start_, zone->wp_, zone);
    if (!s.ok()) return s;
  } else {
    /* For non-sparse files, the data is contigous and we can recover directly
       any missing data using the WP */
    zone->used_capacity_ += to_recover;
    extents_.push_back(new ZoneExtent(extent_start_, to_recover, zone));
  }

  /* Mark up the file as having no missing extents */
  extent_start_ = NO_EXTENT;

  /* Recalculate file size */
  file_size_ = 0;
  for (uint32_t i = 0; i < extents_.size(); i++) {
    file_size_ += extents_[i]->length_;
  }

  return IOStatus::OK();
}

void ZoneFile::ReplaceExtentList(std::vector<ZoneExtent*> new_list) {
  assert(!IsOpenForWR() && new_list.size() > 0);
  assert(new_list.size() == extents_.size());

  WriteLock lck(this);
  extents_ = new_list;
}

void ZoneFile::AddLinkName(const std::string& linkf) {
  linkfiles_.push_back(linkf);
}

IOStatus ZoneFile::RenameLink(const std::string& src, const std::string& dest) {
  auto itr = std::find(linkfiles_.begin(), linkfiles_.end(), src);
  if (itr != linkfiles_.end()) {
    linkfiles_.erase(itr);
    linkfiles_.push_back(dest);
  } else {
    return IOStatus::IOError("RenameLink: Failed to find the linked file");
  }
  return IOStatus::OK();
}

IOStatus ZoneFile::RemoveLinkName(const std::string& linkf) {
  assert(GetNrLinks());
  auto itr = std::find(linkfiles_.begin(), linkfiles_.end(), linkf);
  if (itr != linkfiles_.end()) {
    linkfiles_.erase(itr);
  } else {
    return IOStatus::IOError("RemoveLinkInfo: Failed to find the link file");
  }
  return IOStatus::OK();
}

IOStatus ZoneFile::SetWriteLifeTimeHint(Env::WriteLifeTimeHint lifetime) {
  lifetime_ = lifetime;
  return IOStatus::OK();
}

void ZoneFile::ReleaseActiveZone() {
  assert(active_zone_ != nullptr);
  bool ok = active_zone_->Release();
  assert(ok);
  (void)ok;
  active_zone_ = nullptr;
}

void ZoneFile::SetActiveZone(Zone* zone) {
  assert(active_zone_ == nullptr);
  assert(zone->IsBusy());
  active_zone_ = zone;
}

ZonedWritableFile::ZonedWritableFile(ZonedBlockDevice* zbd, bool _buffered,
                                     std::shared_ptr<ZoneFile> zoneFile,
                                     MetadataWriter* metadata_writer) {
  wp = zoneFile->GetFileSize();

  buffered = _buffered;
  block_sz = zbd->GetBlockSize();
  zoneFile_ = zoneFile;
  buffer_pos = 0;
  sparse_buffer = nullptr;
  buffer = nullptr;

  if (buffered) {
    if (zoneFile->IsSparse()) {
      size_t sparse_buffer_sz;

      sparse_buffer_sz =
          1024 * 1024 + block_sz; /* one extra block size for padding */
      int ret = posix_memalign((void**)&sparse_buffer, sysconf(_SC_PAGESIZE),
                               sparse_buffer_sz);

      if (ret) sparse_buffer = nullptr;

      assert(sparse_buffer != nullptr);

      buffer_sz = sparse_buffer_sz - ZoneFile::SPARSE_HEADER_SIZE - block_sz;
      buffer = sparse_buffer + ZoneFile::SPARSE_HEADER_SIZE;
    } else {
      buffer_sz = 1024 * 1024;
      int ret =
          posix_memalign((void**)&buffer, sysconf(_SC_PAGESIZE), buffer_sz);

      if (ret) buffer = nullptr;
      assert(buffer != nullptr);
    }
  }

  zoneFile_->OpenWR(metadata_writer);
}

ZonedWritableFile::~ZonedWritableFile() {
  IOStatus s = CloseInternal();
  if (buffered) {
    if (sparse_buffer != nullptr) {
      free(sparse_buffer);
    } else {
      free(buffer);
    }
  }

  if (!s.ok()) {
    zoneFile_->GetZbd()->SetZoneDeferredStatus(s);
  }
}

MetadataWriter::~MetadataWriter() {}

IOStatus ZonedWritableFile::Truncate(uint64_t size,
                                     const IOOptions& /*options*/,
                                     IODebugContext* /*dbg*/) {
  zoneFile_->SetFileSize(size);
  return IOStatus::OK();
}

IOStatus ZonedWritableFile::DataSync() {
  if (buffered) {
    IOStatus s;
    buffer_mtx_.lock();
    /* Flushing the buffer will result in a new extent added to the list*/
    s = FlushBuffer();
    buffer_mtx_.unlock();
    if (!s.ok()) {
      return s;
    }

    /* We need to persist the new extent, if the file is not sparse,
     * as we can't use the active zone WP, which is block-aligned, to recover
     * the file size */
    if (!zoneFile_->IsSparse()) return zoneFile_->PersistMetadata();
  } else {
    /* For direct writes, there is no buffer to flush, we just need to push
       an extent for the latest written data */
    zoneFile_->PushExtent();
  }

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::Fsync(const IOOptions& /*options*/,
                                  IODebugContext* /*dbg*/) {
  IOStatus s;
  ZenFSMetricsLatencyGuard guard(zoneFile_->GetZBDMetrics(),
                                 zoneFile_->GetIOType() == IOType::kWAL
                                     ? ZENFS_WAL_SYNC_LATENCY
                                     : ZENFS_NON_WAL_SYNC_LATENCY,
                                 Env::Default());
  zoneFile_->GetZBDMetrics()->ReportQPS(ZENFS_SYNC_QPS, 1);

  s = DataSync();
  if (!s.ok()) return s;

  /* As we've already synced the metadata in DataSync, no need to do it again */
  if (buffered && !zoneFile_->IsSparse()) return IOStatus::OK();

  return zoneFile_->PersistMetadata();
}

IOStatus ZonedWritableFile::Sync(const IOOptions& /*options*/,
                                 IODebugContext* /*dbg*/) {
  return DataSync();
}

IOStatus ZonedWritableFile::Flush(const IOOptions& /*options*/,
                                  IODebugContext* /*dbg*/) {
  return IOStatus::OK();
}

IOStatus ZonedWritableFile::RangeSync(uint64_t offset, uint64_t nbytes,
                                      const IOOptions& /*options*/,
                                      IODebugContext* /*dbg*/) {
  if (wp < offset + nbytes) return DataSync();

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::Close(const IOOptions& /*options*/,
                                  IODebugContext* /*dbg*/) {
  return CloseInternal();
}

IOStatus ZonedWritableFile::CloseInternal() {
  if (!zoneFile_->IsOpenForWR()) {
    return IOStatus::OK();
  }

  IOStatus s = DataSync();
  if (!s.ok()) return s;

  return zoneFile_->CloseWR();
}

IOStatus ZonedWritableFile::FlushBuffer() {
  IOStatus s;

  if (buffer_pos == 0) return IOStatus::OK();

  if (zoneFile_->IsSparse()) {
    s = zoneFile_->SparseAppend(sparse_buffer, buffer_pos);
  } else {
    s = zoneFile_->BufferedAppend(buffer, buffer_pos);
  }

  if (!s.ok()) {
    return s;
  }

  wp += buffer_pos;
  buffer_pos = 0;

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::BufferedWrite(const Slice& slice) {
  uint32_t data_left = slice.size();
  char* data = (char*)slice.data();
  IOStatus s;

  while (data_left) {
    uint32_t buffer_left = buffer_sz - buffer_pos;
    uint32_t to_buffer;

    if (!buffer_left) {
      s = FlushBuffer();
      if (!s.ok()) return s;
      buffer_left = buffer_sz;
    }

    to_buffer = data_left;
    if (to_buffer > buffer_left) {
      to_buffer = buffer_left;
    }

    memcpy(buffer + buffer_pos, data, to_buffer);
    buffer_pos += to_buffer;
    data_left -= to_buffer;
    data += to_buffer;
  }

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::Append(const Slice& data,
                                   const IOOptions& ops/*options*/,
                                   IODebugContext* /*dbg*/) {
  IOStatus s;
  ZenFSMetricsLatencyGuard guard(zoneFile_->GetZBDMetrics(),
                                 zoneFile_->GetIOType() == IOType::kWAL
                                     ? ZENFS_WAL_WRITE_LATENCY
                                     : ZENFS_NON_WAL_WRITE_LATENCY,
                                 Env::Default());
  zoneFile_->GetZBDMetrics()->ReportQPS(ZENFS_WRITE_QPS, 1);
  zoneFile_->GetZBDMetrics()->ReportThroughput(ZENFS_WRITE_THROUGHPUT,
                                               data.size());
  // printf("in ZonedWritableFile::Append, size %ld\n", data.size());
  if (buffered) {
    buffer_mtx_.lock();
    s = BufferedWrite(data);
    buffer_mtx_.unlock();
  } else {
    assert(!ops.metadata);
    s = zoneFile_->Append((void*)data.data(), data.size(), ops.metadata, ops.metadata_num);
    if (s.ok()) wp += data.size();
  }

  return s;
}

IOStatus ZonedWritableFile::PositionedAppend(const Slice& data, uint64_t offset,
                                             const IOOptions& ops/*options*/,
                                             IODebugContext* /*dbg*/) {
  IOStatus s;
  ZenFSMetricsLatencyGuard guard(zoneFile_->GetZBDMetrics(),
                                 zoneFile_->GetIOType() == IOType::kWAL
                                     ? ZENFS_WAL_WRITE_LATENCY
                                     : ZENFS_NON_WAL_WRITE_LATENCY,
                                 Env::Default());
  zoneFile_->GetZBDMetrics()->ReportQPS(ZENFS_WRITE_QPS, 1);
  zoneFile_->GetZBDMetrics()->ReportThroughput(ZENFS_WRITE_THROUGHPUT,
                                               data.size());
#ifndef NDBG
  printf("in ZonedWritableFile::PositionedAppend, offset %ld wp %ld size %ld\n", offset, wp, data.size());
#endif
  if (offset != wp) {
    return IOStatus::IOError("positioned append not at write pointer");
  }

  if (buffered) {
    assert(false);
    buffer_mtx_.lock();
    s = BufferedWrite(data);
    buffer_mtx_.unlock();
  } else {
    size_t asize = data.size();
    s = zoneFile_->Append((void*)data.data(), data.size(), ops.metadata, ops.metadata_num, &asize);
    if (s.ok()) wp += asize;
  }

  return s;
}

void ZonedWritableFile::SetWriteLifeTimeHint(Env::WriteLifeTimeHint hint) {
  zoneFile_->SetWriteLifeTimeHint(hint);
}

IOStatus ZonedSequentialFile::Read(size_t n, const IOOptions& opts /*options*/,
                                   Slice* result, char* scratch,
                                   IODebugContext* /*dbg*/) {
  IOStatus s;

  s = zoneFile_->PositionedRead(rp, n, result, scratch, opts);
  if (s.ok()) rp += result->size();

  return s;
}

IOStatus ZonedSequentialFile::Skip(uint64_t n) {
  if (rp + n >= zoneFile_->GetFileSize())
    return IOStatus::InvalidArgument("Skip beyond end of file");
  rp += n;
  return IOStatus::OK();
}

IOStatus ZonedSequentialFile::PositionedRead(uint64_t offset, size_t n,
                                             const IOOptions& opts /*options*/,
                                             Slice* result, char* scratch,
                                             IODebugContext* /*dbg*/) {
  assert(false);
  return zoneFile_->PositionedRead(offset, n, result, scratch, opts);
}

IOStatus ZonedRandomAccessFile::Read(uint64_t offset, size_t n,
                                     const IOOptions& opts /*options*/,
                                     Slice* result, char* scratch,
                                     IODebugContext* /*dbg*/) const {
  return zoneFile_->PositionedRead(offset, n, result, scratch, opts);
}

size_t ZoneFile::GetUniqueId(char* id, size_t max_size) {
  /* Based on the posix fs implementation */
  if (max_size < kMaxVarint64Length * 3) {
    return 0;
  }

  struct stat buf;
  int fd = zbd_->GetReadFD();
  int result = fstat(fd, &buf);
  if (result == -1) {
    return 0;
  }

  char* rid = id;
  rid = EncodeVarint64(rid, buf.st_dev);
  rid = EncodeVarint64(rid, buf.st_ino);
  rid = EncodeVarint64(rid, file_id_);
  assert(rid >= id);
  return static_cast<size_t>(rid - id);

  return 0;
}

IOStatus ZoneFile::MigrateData(uint64_t offset, uint32_t length,
                               Zone* target_zone) {
  uint32_t step = 128 << 10;
  uint32_t read_sz = step;
  int block_sz = zbd_->GetBlockSize();
  printf("in ZoneFile::MigrateData\n");
  assert(offset % block_sz == 0);
  if (offset % block_sz != 0) {
    return IOStatus::IOError("MigrateData offset is not aligned!\n");
  }

  char* buf;
  int ret = posix_memalign((void**)&buf, block_sz, step);
  if (ret) {
    return IOStatus::IOError("failed allocating alignment write buffer\n");
  }

  int pad_sz = 0;
  while (length > 0) {
    read_sz = length > read_sz ? read_sz : length;
    pad_sz = read_sz % block_sz == 0 ? 0 : (block_sz - (read_sz % block_sz));

    int r = zbd_->DirectRead(buf, offset, read_sz + pad_sz);
    if (r < 0) {
      free(buf);
      return IOStatus::IOError(strerror(errno));
    }
    target_zone->Append(buf, r);
    length -= read_sz;
    offset += r;
  }

  free(buf);

  return IOStatus::OK();
}

size_t ZonedRandomAccessFile::GetUniqueId(char* id, size_t max_size) const {
  return zoneFile_->GetUniqueId(id, max_size);
}

}  // namespace ROCKSDB_NAMESPACE

#endif  // !defined(ROCKSDB_LITE) && !defined(OS_WIN)
