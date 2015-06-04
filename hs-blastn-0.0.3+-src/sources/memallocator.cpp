#include "memallocator.h"

void Chunk::Clear(std::size_t blockSize, Uint1 blocks)
{
    firstAvailableBlock_ = 0;
    blocksAvailable_ = blocks;
    Uint1 i = 0;
    Uint1* p = pData_;
    for (; i != blocks; p += blockSize)
    {
        *p = ++i;
    }
}

void Chunk::Init(std::size_t blockSize, Uint1 blocks)
{
    pData_ = new Uint1[blockSize * blocks];
    firstAvailableBlock_ = 0;
    blocksAvailable_ = blocks;
    Uint1 i = 0;
    Uint1* p = pData_;
    for (; i != blocks; p += blockSize)
    {
        *p = ++i;
    }
}

void* Chunk::Allocate(std::size_t blockSize)
{
    if (!blocksAvailable_) return NULL;
    Uint1* pResult =
        pData_ + (firstAvailableBlock_ * blockSize);
    firstAvailableBlock_ = *pResult;
    --blocksAvailable_;
    return pResult;
}

void Chunk::Deallocate(void* p, std::size_t blockSize)
{
    assert(p >= pData_);
    Uint1* toRelease = static_cast<Uint1*>(p);
    // Alignment check
    assert((toRelease - pData_) % blockSize == 0);
    *toRelease = firstAvailableBlock_;
    firstAvailableBlock_ = static_cast<Uint1>((toRelease - pData_) / blockSize);
    // Truncation check
    assert(firstAvailableBlock_ ==
        (toRelease -pData_) / blockSize);
    ++blocksAvailable_;
}

void Chunk::Release()
{
    if (pData_)
    {
        delete[] pData_;
        firstAvailableBlock_ = blocksAvailable_ = 0;
        pData_ = NULL;
    }
}

void FixedAllocator::Release()
{
    std::vector<Chunk>::iterator iter;
    for (iter = chunks_.begin(); iter != chunks_.end(); ++iter)
    {
        iter->Release();
    }
}

void FixedAllocator::Init(std::size_t blockSize, Uint1 blocks)
{
    blockSize_ = blockSize;
    numBlocks_ = blocks;
    curChunk_ = 0;

    Chunk new_chunk;
    new_chunk.Init(blockSize, blocks);
    chunks_.push_back(new_chunk);

    deallocIndex_ = 0;
}

void FixedAllocator::Clear()
{
    curChunk_ = 0;
    deallocIndex_ = 0;
    chunks_[curChunk_].Clear(blockSize_, numBlocks_);
}

void FixedAllocator::GetOneChunk()
{
    if (curChunk_ + 1< chunks_.size())
    {
        ++curChunk_;
        chunks_[curChunk_].Clear(blockSize_, numBlocks_);
    }
    else
    {
        Chunk chunk;
        chunk.Init(blockSize_, numBlocks_);
        chunks_.push_back(chunk);
        ++curChunk_;
    }
}

void* FixedAllocator::Allocate()
{
    if (deallocIndex_ > 0)
        return deallocate_[--deallocIndex_];

    if (!chunks_[curChunk_].blocksAvailable_)
    {
        GetOneChunk();
    }

    return chunks_[curChunk_].Allocate(blockSize_);
}

void FixedAllocator::Deallocate(void* p)
{
    if (deallocIndex_ < kMaxDeallocate)
        deallocate_[deallocIndex_++] = p;
}


void SmallObjAllocator::Release()
{
    std::vector<FixedAllocator>::iterator iter;
    for (iter = pool_.begin(); iter != pool_.end(); ++iter)
    {
        iter->Release();
    }
}

void*
SmallObjAllocator::Allocate(std::size_t numBytes)
{
    if (numBytes > kSoaMaxObjectSize)
    {
        unsigned char* pResult = new unsigned char[numBytes];
        return pResult;
    }

    std::size_t idx = (numBytes - 1) / kSoaIncObjectSize;
    return pool_[idx].Allocate();
}

void
SmallObjAllocator::Deallocate(void* p, std::size_t size)
{
    if (size > kSoaMaxObjectSize)
    {
        unsigned char* q = static_cast<unsigned char*>(p);
        delete[] q;
        return;
    }

    std::size_t idx = (size - 1) / kSoaIncObjectSize;
    pool_[idx].Deallocate(p);
}

void SmallObjAllocator::Clear()
{
    std::vector<FixedAllocator>::iterator iter;
    for (iter = pool_.begin(); iter != pool_.end(); ++iter)
    {
        iter->Clear();
    }
}

SmallObjAllocator::SmallObjAllocator()
{
    std::size_t size = kSoaMinObjectSize;
    std::size_t n;
    for (; size <= kSoaMaxObjectSize; size += kSoaIncObjectSize)
    {
        n = kSoaChunkSize / size;
        assert (n <= UINT1_MAX);
        FixedAllocator fa;
        fa.Init(size, n);
        pool_.push_back(fa);
    }
}
