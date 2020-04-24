#ifndef BUFFER_LINE_ITERATOR_H
#define BUFFER_LINE_ITERATOR_H

#include <fstream>
#include <iterator>
#include <sstream>
#include <string>

#include "def.h"

class BufferLineIterator
{
public:
    typedef std::input_iterator_tag     iterator_tag;
    typedef std::string                 value_type;
    typedef std::ptrdiff_t              difference_type;
    typedef const std::string*               pointer;
    typedef const std::string&          reference;

public:
    BufferLineIterator() : ins_(NULL), is_valid_(false) { buf_ = new char[kBufferSize]; x_InitStreamState(); }

    BufferLineIterator(std::ifstream& s) : ins_(&s), is_valid_( (*ins_) ? true : false) { buf_ = new char[kBufferSize]; x_InitStreamState(); }
	
	~BufferLineIterator() { delete[] buf_; }

    reference operator*() { return line_; }

    pointer operator->() const { return &line_; }

    BufferLineIterator& operator++()
    {
        x_GetOneLine();
        return *this;
    }

    //BufferLineIterator operator++(int)
    //{
    //    BufferLineIterator tmp = *this;
    //    x_GetOneLine();
    //    return tmp;
    //}

    void reset(std::ifstream* new_stream)
    {
        ins_ = new_stream;
        if (*ins_)
            is_valid_ = true;
        x_InitStreamState();
    }

    bool eof()
    {
        return !is_valid_;
    }

    void clear()
    {
        ins_ = NULL;
        line_.clear();
        is_valid_ = false;
        x_InitStreamState();
    }

    std::string& get_string()
    {
        return line_;
    }

private:
    bool x_GetOneLine();
    bool x_ReadBuffer();
    void x_InitStreamState();

private:
    std::ifstream*           ins_;
    static const index_t    kBufferSize = 1024 * 1024;
    char*                   buf_;
    index_t                 cur_;
    index_t                 buf_sz_;
    bool                    done_;

private:
    std::string line_;
    bool        is_valid_;
};

#endif //BUFFER_LINE_ITERATOR_H
