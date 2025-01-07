//  Copyright (c) 2011-present, Facebook, Inc.  All rights reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
//
// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "table/block_fetcher.h"

#include <cassert>
#include <cinttypes>
#include <string>

#include "logging/logging.h"
#include "memory/memory_allocator_impl.h"
#include "monitoring/perf_context_imp.h"
#include "rocksdb/compression_type.h"
#include "rocksdb/env.h"
#include "table/block_based/block.h"
#include "table/block_based/block_based_table_reader.h"
#include "table/block_based/block_type.h"
#include "table/block_based/reader_common.h"
#include "table/format.h"
#include "table/persistent_cache_helper.h"
#include "util/compression.h"
#include "util/stop_watch.h"

namespace ROCKSDB_NAMESPACE {

inline void BlockFetcher::ProcessTrailerIfPresent() {
  if (footer_.GetBlockTrailerSize() > 0) {
    assert(footer_.GetBlockTrailerSize() == BlockBasedTable::kBlockTrailerSize);
    if (read_options_.verify_checksums) {
      io_status_ = status_to_io_status(
          VerifyBlockChecksum(footer_, slice_.data(), block_size_,
                              file_->file_name(), handle_.offset()));
      RecordTick(ioptions_.stats, BLOCK_CHECKSUM_COMPUTE_COUNT);
      if (!io_status_.ok()) {
        assert(io_status_.IsCorruption());
        RecordTick(ioptions_.stats, BLOCK_CHECKSUM_MISMATCH_COUNT);
      }
    }
    compression_type_ =
        BlockBasedTable::GetBlockCompressionType(slice_.data(), block_size_);
  } else {
    // E.g. plain table or cuckoo table
    compression_type_ = kNoCompression;
  }
}

inline void BlockFetcher::ProcessTrailerIfPresent_CZNS(size_t asize) {
  assert(footer_.GetBlockTrailerSize() == BlockBasedTable::kBlockTrailerSize);

  size_t block_size = 0;
  if (prefetch_buffer_ != nullptr) {
    // Use hash_map to fetch the correct decompressed data
    auto& hash_map = prefetch_buffer_->GetFirstBuffer()->hash_map; // TO FIX
    auto it = hash_map.find(handle_.offset());
    if (it != hash_map.end()) {
      assert(!asize);
#ifndef NDBG
      printf("in ProcessTrailerIfPresent_CZNS, bs %d\n", it->second.second);
#endif
      block_size = it->second.second;
    }
  } else {
    assert(asize);
#ifndef NDBG
    printf("in ProcessTrailerIfPresent_CZNS, bs %lu\n", asize);
#endif
    block_size = asize;
  }
  block_size -= footer_.GetBlockTrailerSize();
  compression_type_ =
    BlockBasedTable::GetBlockCompressionType(slice_.data(), block_size);
}

inline bool BlockFetcher::TryGetUncompressBlockFromPersistentCache() {
  if (cache_options_.persistent_cache &&
      !cache_options_.persistent_cache->IsCompressed()) {
    Status status = PersistentCacheHelper::LookupUncompressed(
        cache_options_, handle_, contents_);
    if (status.ok()) {
      // uncompressed page is found for the block handle
      return true;
    } else {
      // uncompressed page is not found
      if (ioptions_.logger && !status.IsNotFound()) {
        assert(!status.ok());
        ROCKS_LOG_INFO(ioptions_.logger,
                       "Error reading from persistent cache. %s",
                       status.ToString().c_str());
      }
    }
  }
  return false;
}

inline bool BlockFetcher::TryGetFromPrefetchBuffer() {
#ifndef NDBG
  printf("in BlockFetcher::TryGetFromPrefetchBuffer\n");
#endif
  if (prefetch_buffer_ != nullptr) {
    IOOptions opts;
    IOStatus io_s = file_->PrepareIOOptions(read_options_, opts);
    if (io_s.ok()) {
      // RCZNS compresses data along with trailer and decompressed them together, 
      // thus here we do NOT need to add block size with trailer size.
      uint64_t alens = 0;
      opts.alens_addr = &alens;
      opts.decompressed_read = decompressed_read_;
      opts.flexible_mode = true;
      opts.actual_size = handle_.prev_size();

      bool read_from_prefetch_buffer = prefetch_buffer_->TryReadFromCache(
          opts, file_, handle_.offset(), block_size_, &slice_,
          &io_s, for_compaction_);
      if (read_from_prefetch_buffer) {
#ifndef NDBG
        printf("read_from_prefetch_buffer, continue\n");
#endif
        if (opts.decompressed_read)
          ProcessTrailerIfPresent_CZNS();
        else
          ProcessTrailerIfPresent();
        if (!io_status_.ok()) {
          return true;
        }
        got_from_prefetch_buffer_ = true;
        used_buf_ = const_cast<char*>(slice_.data());
#ifndef NDBG
        printf("used_bufsz %lu\n", slice_.size());
#endif
      }
    }
    if (!io_s.ok()) {
      io_status_ = io_s;
      return true;
    }
  }
  return got_from_prefetch_buffer_;
}

inline bool BlockFetcher::TryGetSerializedBlockFromPersistentCache() {
  if (cache_options_.persistent_cache &&
      cache_options_.persistent_cache->IsCompressed()) {
    assert(false);
    std::unique_ptr<char[]> buf;
    io_status_ = status_to_io_status(PersistentCacheHelper::LookupSerialized(
        cache_options_, handle_, &buf, block_size_with_trailer_));
    if (io_status_.ok()) {
      heap_buf_ = CacheAllocationPtr(buf.release());
      used_buf_ = heap_buf_.get();
      slice_ = Slice(heap_buf_.get(), block_size_);
      ProcessTrailerIfPresent();
      return true;
    } else if (!io_status_.IsNotFound() && ioptions_.logger) {
      assert(!io_status_.ok());
      ROCKS_LOG_INFO(ioptions_.logger,
                     "Error reading from persistent cache. %s",
                     io_status_.ToString().c_str());
    }
  }
  return false;
}

inline void BlockFetcher::PrepareBufferForBlockFromFile() {
  // cache miss read from device
  if ((do_uncompress_ || ioptions_.allow_mmap_reads) &&
      block_size_with_trailer_ < kDefaultStackBufferSize) {
    // If we've got a small enough chunk of data, read it in to the
    // trivially allocated stack buffer instead of needing a full malloc()
    //
    // `GetBlockContents()` cannot return this data as its lifetime is tied to
    // this `BlockFetcher`'s lifetime. That is fine because this is only used
    // in cases where we do not expect the `GetBlockContents()` result to be the
    // same buffer we are assigning here. If we guess incorrectly, there will be
    // a heap allocation and memcpy in `GetBlockContents()` to obtain the final
    // result. Considering we are eliding a heap allocation here by using the
    // stack buffer, the cost of guessing incorrectly here is one extra memcpy.
    //
    // When `do_uncompress_` is true, we expect the uncompression step will
    // allocate heap memory for the final result. However this expectation will
    // be wrong if the block turns out to already be uncompressed, which we
    // won't know for sure until after reading it.
    //
    // When `ioptions_.allow_mmap_reads` is true, we do not expect the file
    // reader to use the scratch buffer at all, but instead return a pointer
    // into the mapped memory. This expectation will be wrong when using a
    // file reader that does not implement mmap reads properly.
    used_buf_ = &stack_buf_[0];
  } else if (maybe_compressed_ && !do_uncompress_) {
    compressed_buf_ =
        AllocateBlock(block_size_with_trailer_, memory_allocator_compressed_);
    used_buf_ = compressed_buf_.get();
  } else {
    heap_buf_ = AllocateBlock(block_size_with_trailer_, memory_allocator_);
    used_buf_ = heap_buf_.get();
  }
}

inline void BlockFetcher::InsertCompressedBlockToPersistentCacheIfNeeded() {
  if (io_status_.ok() && read_options_.fill_cache &&
      cache_options_.persistent_cache &&
      cache_options_.persistent_cache->IsCompressed()) {
    PersistentCacheHelper::InsertSerialized(cache_options_, handle_, used_buf_,
                                            block_size_with_trailer_);
  }
}

inline void BlockFetcher::InsertUncompressedBlockToPersistentCacheIfNeeded() {
  if (io_status_.ok() && !got_from_prefetch_buffer_ &&
      read_options_.fill_cache && cache_options_.persistent_cache &&
      !cache_options_.persistent_cache->IsCompressed()) {
    // insert to uncompressed cache
    PersistentCacheHelper::InsertUncompressed(cache_options_, handle_,
                                              *contents_);
  }
}

inline void BlockFetcher::CopyBufferToHeapBuf(size_t asize) {
#ifndef NDBG
  printf("in BlockFetcher::CopyBufferToHeapBuf, dr? %d\n", decompressed_read_);
#endif
  size_t sz = decompressed_read_? asize : block_size_with_trailer_;
  assert(used_buf_ != heap_buf_.get());
  heap_buf_ = AllocateBlock(sz, memory_allocator_);
  memcpy(heap_buf_.get(), used_buf_, sz);
#ifndef NDEBUG
  num_heap_buf_memcpy_++;
#endif
}

inline void BlockFetcher::CopyBufferToCompressedBuf() {
  assert(used_buf_ != compressed_buf_.get());
  compressed_buf_ =
      AllocateBlock(block_size_with_trailer_, memory_allocator_compressed_);
  memcpy(compressed_buf_.get(), used_buf_, block_size_with_trailer_);
#ifndef NDEBUG
  num_compressed_buf_memcpy_++;
#endif
}

// Entering this method means the block is not compressed or do not need to be
// uncompressed. The block can be in one of the following buffers:
// 1. prefetch buffer if prefetch is enabled and the block is prefetched before
// 2. stack_buf_ if block size is smaller than the stack_buf_ size and block
//    is not compressed
// 3. heap_buf_ if the block is not compressed
// 4. compressed_buf_ if the block is compressed
// 5. direct_io_buf_ if direct IO is enabled
// After this method, if the block is compressed, it should be in
// compressed_buf_, otherwise should be in heap_buf_.
inline void BlockFetcher::GetBlockContents() {
#ifndef NDBG
  printf("in BlockFetcher::GetBlockContents slice_sz %lu bswt %lu dr? %d\n", slice_.size(), 
          block_size_with_trailer_, decompressed_read_);
#endif
  size_t content_size = decompressed_read_? (slice_.size() - footer_.GetBlockTrailerSize()) : block_size_;
  if (slice_.data() != used_buf_) {
    // the slice content is not the buffer provided
    *contents_ = BlockContents(Slice(slice_.data(), content_size));
  } else {
    // page can be either uncompressed or compressed, the buffer either stack
    // or heap provided. Refer to https://github.com/facebook/rocksdb/pull/4096
    if (got_from_prefetch_buffer_ || used_buf_ == &stack_buf_[0]) {
      CopyBufferToHeapBuf(slice_.size());
    } else if (used_buf_ == compressed_buf_.get()) {
      assert(false);
      if (compression_type_ == kNoCompression &&
          memory_allocator_ != memory_allocator_compressed_) {
        CopyBufferToHeapBuf();
      } else {
        heap_buf_ = std::move(compressed_buf_);
      }
    } else if (direct_io_buf_.get() != nullptr) {
      if (compression_type_ == kNoCompression) {
        CopyBufferToHeapBuf(slice_.size());
      } else {
        assert(false);
        CopyBufferToCompressedBuf();
        heap_buf_ = std::move(compressed_buf_);
      }
    }
    *contents_ = BlockContents(std::move(heap_buf_), content_size);
  }
#ifndef NDEBUG
  contents_->has_trailer = footer_.GetBlockTrailerSize() > 0;
#endif
}

IOStatus BlockFetcher::ReadBlockContents() {
#ifndef NDBG
  printf("in BlockFetcher::ReadBlockContents, direct? %d dr? %d fx? %d\n", file_->use_direct_io(), 
          decompressed_read_, flexible_mode_);
#endif
  if (TryGetUncompressBlockFromPersistentCache()) {
    compression_type_ = kNoCompression;
#ifndef NDEBUG
    contents_->has_trailer = footer_.GetBlockTrailerSize() > 0;
#endif  // NDEBUG
    return IOStatus::OK();
  }
  if (TryGetFromPrefetchBuffer()) {
    if (!io_status_.ok()) {
      return io_status_;
    }
  } else if (!TryGetSerializedBlockFromPersistentCache()) {
    IOOptions opts;
    uint64_t alens = 0;
    size_t asize = block_size_with_trailer_;
    if (decompressed_read_) {
      opts.alens_addr = &alens;
      opts.decompressed_read = decompressed_read_;
      opts.flexible_mode = flexible_mode_;
      opts.actual_size = handle_.prev_size();
      asize = block_size_;
    } 
    io_status_ = file_->PrepareIOOptions(read_options_, opts);
    // Actual file read
    if (io_status_.ok()) {
      if (file_->use_direct_io()) {
        PERF_TIMER_GUARD(block_read_time);
        PERF_CPU_TIMER_GUARD(
            block_read_cpu_time,
            ioptions_.env ? ioptions_.env->GetSystemClock().get() : nullptr);
        io_status_ =
            file_->Read(opts, handle_.offset(), asize, // block_size_with_trailer_,
                        &slice_, nullptr, &direct_io_buf_);
        PERF_COUNTER_ADD(block_read_count, 1);
        used_buf_ = const_cast<char*>(slice_.data());
        asize = slice_.size();
#ifndef NDBG
        printf("received buf size %lu\n", slice_.size());
#endif
      } else {
        assert(false);
        PrepareBufferForBlockFromFile();
        PERF_TIMER_GUARD(block_read_time);
        PERF_CPU_TIMER_GUARD(
            block_read_cpu_time,
            ioptions_.env ? ioptions_.env->GetSystemClock().get() : nullptr);
        io_status_ =
            file_->Read(opts, handle_.offset(), block_size_with_trailer_,
                        &slice_, used_buf_, nullptr);
        PERF_COUNTER_ADD(block_read_count, 1);
#ifndef NDEBUG
        if (slice_.data() == &stack_buf_[0]) {
          num_stack_buf_memcpy_++;
        } else if (slice_.data() == heap_buf_.get()) {
          num_heap_buf_memcpy_++;
        } else if (slice_.data() == compressed_buf_.get()) {
          num_compressed_buf_memcpy_++;
        }
#endif
      }
    }
    
    // The after-compressed data length is not the same as the after-DEcompressed
    // data length, thus we need a hash_map to record the correct locations in 
    // the prefetch buffer.
    if (prefetch_buffer_ != nullptr) {
      uint16_t* tmp = nullptr;
      auto buf = prefetch_buffer_->GetFirstBuffer();
      auto& hash_map = buf->hash_map;

      if (opts.alens_addr) {
        assert(opts.decompressed_read);
        tmp = (uint16_t *)(*opts.alens_addr);
        hash_map.clear();
        uint64_t offset_ac = handle_.offset(), offset_ad = handle_.offset() - buf->offset_;
#ifndef NDBG
        uint64_t min_offset_ac = offset_ac, min_offset_ad = offset_ad;
#endif
        for (int i = 1; i <= *tmp; i++) {
          offset_ac += tmp[i << 1];
          offset_ad += tmp[(i << 1) - 1];
        }
#ifndef NDBG
        printf("renew hashmap before decompressed [%lu, %lu] -> after [%lu, %lu], cursz %lu ofst %lu\n", min_offset_ac, offset_ac, min_offset_ad, offset_ad,
                hash_map.size(), buf->offset_);
#endif
        free(tmp);
      }
    }

    // TODO: introduce dedicated perf counter for range tombstones
    switch (block_type_) {
      case BlockType::kFilter:
      case BlockType::kFilterPartitionIndex:
        PERF_COUNTER_ADD(filter_block_read_count, 1);
        break;

      case BlockType::kCompressionDictionary:
        PERF_COUNTER_ADD(compression_dict_block_read_count, 1);
        break;

      case BlockType::kIndex:
        PERF_COUNTER_ADD(index_block_read_count, 1);
        break;

      // Nothing to do here as we don't have counters for the other types.
      default:
        break;
    }

    PERF_COUNTER_ADD(block_read_byte, block_size_with_trailer_);
    if (!io_status_.ok()) {
      return io_status_;
    }

    if (!decompressed_read_ && slice_.size() != block_size_with_trailer_) {
      return IOStatus::Corruption(
          "truncated block read from " + file_->file_name() + " offset " +
          std::to_string(handle_.offset()) + ", expected " +
          std::to_string(block_size_with_trailer_) + " bytes, got " +
          std::to_string(slice_.size()));
    }
    if (decompressed_read_)
      ProcessTrailerIfPresent_CZNS(asize);
    else
      ProcessTrailerIfPresent();
    if (io_status_.ok()) {
      InsertCompressedBlockToPersistentCacheIfNeeded();
    } else {
      return io_status_;
    }
  }

  if (do_uncompress_ && compression_type_ != kNoCompression) {
    PERF_TIMER_GUARD(block_decompress_time);
    // compressed page, uncompress, update cache
    UncompressionContext context(compression_type_);
    UncompressionInfo info(context, uncompression_dict_, compression_type_);
    io_status_ = status_to_io_status(UncompressSerializedBlock(
        info, slice_.data(), block_size_, contents_, footer_.format_version(),
        ioptions_, memory_allocator_));
#ifndef NDEBUG
    num_heap_buf_memcpy_++;
#endif
    // Save the compressed block without trailer
    slice_ = Slice(slice_.data(), block_size_);
  } else {
    GetBlockContents();
    slice_ = Slice();
  }

  InsertUncompressedBlockToPersistentCacheIfNeeded();

  return io_status_;
}

IOStatus BlockFetcher::ReadAsyncBlockContents() {
  assert(false);
  if (TryGetUncompressBlockFromPersistentCache()) {
    compression_type_ = kNoCompression;
#ifndef NDEBUG
    contents_->has_trailer = footer_.GetBlockTrailerSize() > 0;
#endif  // NDEBUG
    return IOStatus::OK();
  } else if (!TryGetSerializedBlockFromPersistentCache()) {
    assert(prefetch_buffer_ != nullptr);
    if (!for_compaction_) {
      IOOptions opts;
      IOStatus io_s = file_->PrepareIOOptions(read_options_, opts);
      if (!io_s.ok()) {
        return io_s;
      }
      io_s = status_to_io_status(prefetch_buffer_->PrefetchAsync(
          opts, file_, handle_.offset(), block_size_with_trailer_, &slice_));
      if (io_s.IsTryAgain()) {
        return io_s;
      }
      if (io_s.ok()) {
        // Data Block is already in prefetch.
        got_from_prefetch_buffer_ = true;
        ProcessTrailerIfPresent();
        if (!io_status_.ok()) {
          return io_status_;
        }
        used_buf_ = const_cast<char*>(slice_.data());

        if (do_uncompress_ && compression_type_ != kNoCompression) {
          PERF_TIMER_GUARD(block_decompress_time);
          // compressed page, uncompress, update cache
          UncompressionContext context(compression_type_);
          UncompressionInfo info(context, uncompression_dict_,
                                 compression_type_);
          io_status_ = status_to_io_status(UncompressSerializedBlock(
              info, slice_.data(), block_size_, contents_,
              footer_.format_version(), ioptions_, memory_allocator_));
#ifndef NDEBUG
          num_heap_buf_memcpy_++;
#endif
        } else {
          GetBlockContents();
        }
        InsertUncompressedBlockToPersistentCacheIfNeeded();
        return io_status_;
      }
    }
    // Fallback to sequential reading of data blocks in case of io_s returns
    // error or for_compaction_is true.
    return ReadBlockContents();
  }
  return io_status_;
}

}  // namespace ROCKSDB_NAMESPACE
