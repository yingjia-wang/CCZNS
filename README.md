## Reviving In-Storage Hardware Compression on ZNS SSDs through Host-SSD Collaboration (HPCA'25)

This repo contains the source code of our paper, including:
- FEMU that implements CCZNS
- CCZNS-aware Linux kernel
- CCZNS-aware RocksDB and ZenFS

### Notes
FEMU is augmented with computational support to emulate hardware compression and the new I/O processing architecture is introduced in Sec V-A. 
Notably, our modified FEMU necessitates more CPU cores (e.g., 32 in our evaluation, see Sec V-A) and higher memory in the server.
The reason for requiring higher memory is that, to ensure the slow CPU-emulated decompression does not block the fast read I/O, FEMU actually stores the before-compressed data (which is typically larger) in the reserved server memory so that subsequent reads can directly fetch the data without decompression.

