#include "buffer_line_iterator.h"

bool BufferLineIterator::x_GetOneLine()
{
    line_.clear();
    if (eof()) return false;
    while (1)
    {
        const char* start = buf_ + cur_;
        const char* end = buf_ + buf_sz_;
        const char* p = start;
        for (p = start; p < end; ++p)
        {
            if ( *p == '\n')
            {
                cur_ += p - start + 1;
                break;
            }
        }
        line_.append(start, p - start);
        if (p < end) break;
        if (p == end && !x_ReadBuffer()) break;
    }

    if (line_.size() > 0 && line_[line_.size() - 1] == '\r') line_.resize(line_.size() - 1);
    return true;
}


bool BufferLineIterator::x_ReadBuffer()
{
    if (done_) { is_valid_ = false; return false; }
    std::streambuf* sb = ins_->rdbuf();
    const bool ok = sb && ins_->good();
    std::streamsize r = ok ? sb->sgetn(buf_, kBufferSize) : 0;
    buf_sz_ = static_cast<index_t>(r);
    cur_ = 0;
    bool ret = true;
    if (buf_sz_ == 0)
    {
        done_ = true;
        ret = false;
    }
    else if (buf_sz_ < kBufferSize)
    {
        done_ = true;
    }
    if ( !ret ) is_valid_ = false;
    return ret;
}

void BufferLineIterator::x_InitStreamState()
{
    cur_ = buf_sz_ = kBufferSize;
    done_ = false;
}

