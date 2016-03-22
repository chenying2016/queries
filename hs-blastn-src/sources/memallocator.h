/*
 * memory allocator
 * see "Modern C++ Design: Generic Programming and Design Patterns Applied" by Andrei Alexandrescu
*/

#ifndef MEMALLOCATOR_H
#define	MEMALLOCATOR_H

#include <vector>
#include "def.h"

struct Chunk
{
    void Init(std::size_t blockSize, Uint1 blocks);
    void Release();
    void* Allocate(std::size_t blockSize);
    void Deallocate(void* p, std::size_t blockSize);
    void Clear(std::size_t blockSize, Uint1 blocks);

    Uint1* pData_;
    Uint1 firstAvailableBlock_,
          blocksAvailable_;
};

class FixedAllocator
{
private:
    std::size_t blockSize_;
    Uint1 numBlocks_;
    typedef std::vector<Chunk> Chunks;
    Chunks chunks_;
    std::size_t curChunk_;

    static const std::size_t kMaxDeallocate = 128;
    void* deallocate_[kMaxDeallocate];
    std::size_t deallocIndex_;

public:
    void GetOneChunk();
    void* Allocate();
    void Deallocate(void* p);
    void Release();
    void Init(std::size_t blockSize, Uint1 blocks);
    void Clear();
};

class SmallObjAllocator
{
public:
    static const std::size_t kSoaMaxObjectSize = 128;
    static const std::size_t kSoaMinObjectSize = 32;
    static const std::size_t kSoaIncObjectSize = 32;
    static const std::size_t kSoaChunkSize = 4096;
public:
    SmallObjAllocator();
    void* Allocate(std::size_t numBytes);
    void Deallocate(void* p, std::size_t size);
    void Clear();
    void Release();

private:
    std::vector<FixedAllocator> pool_;
};

#endif	/* MEMALLOCATOR_H */

