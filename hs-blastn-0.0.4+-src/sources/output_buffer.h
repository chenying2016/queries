#ifndef RESULTS_BUFFER_H_
#define RESULTS_BUFFER_H_

#include <string>
#include "def.h"
#include "utility.h"

using std::string;

class ResultBlk
{
public:
	void clear()
	{
		_current_position = 0;
		_left_size = kBlkSize;
	}
	void init()
	{
		_space = (Uint1*)calloc(1, kBlkSize);
		_next_block = NULL;
		clear();
	}
	void destroy()
	{
		free(_space);
		_space = NULL;
		clear();
	}
	ResultBlk()
	{
		_space = NULL;
		init();
	}
	~ResultBlk()
	{
		destroy();
	}
	void set_next_block(ResultBlk* next_blk)
	{
		_next_block = next_blk;
	}
	ResultBlk* get_next_block()
	{
		return _next_block;
	}
	
public:
	bool check_space_size(Uint4 size)
	{
		return (_current_position + size) < kBlkSize;
	}
	void mark_end()
	{
		_space[_current_position] = '\0';
	}
	void flush(FILE* file)
	{
		const char* out = (const char*)_space;
		mark_end();
		if (file == NULL)
		{
			fprintf(stdout, "%s", out);
		}
		else
		{
			fprintf(file, "%s", out);
		}
	}
	void append_data(const char* str, Uint4 size)
	{
		memcpy(_space + _current_position, str, size);
		_current_position += size;
	}
	void append_data(const char* str)
	{
		if (str == NULL) return;
		append_data(str);
	}
	void append_data(char ch)
	{
		_space[_current_position] = ch;
		++_current_position;
	}
	
private:
	static const Uint4 kBlkSize = (1ULL<<20);
	Uint1* _space;
	ResultBlk* _next_block;
	Uint4  _current_position;
	Uint4  _left_size;
};

class ResultsBuffer
{
public:
	ResultsBuffer()
	{
		blocks = new ResultBlk;
		current_block = blocks;
	}
	~ResultsBuffer()
	{
		Destroy();
	}
	void clear()
	{
		current_block = blocks;
		current_block->clear();
	}
	void advance_block()
	{
		ResultBlk* next = current_block->get_next_block();
		if (next == NULL)
			next = new ResultBlk;
		current_block->set_next_block(next);
		current_block = next;
		current_block->clear();
	}
	void flush(FILE* file)
	{
		current_block->mark_end();
		ResultBlk* p = blocks;
		while (p != current_block->get_next_block())
		{
			p->flush(file);
			p = p->get_next_block();
		}
	}
	void Destroy()
	{
		ResultBlk* next;
		while (blocks != NULL)
		{
			next = blocks->get_next_block();
			delete blocks;
			blocks = next;
		}
		blocks = NULL;
	}
	void Init(const char*) {}
	void CloseFile() {}
	
public:
	ResultsBuffer& operator << (char ch)
	{
		if (!current_block->check_space_size(sizeof(ch)))
		{
			advance_block();
		}
		current_block->append_data(ch);
		return *this;
	}
	ResultsBuffer& operator << (const char* str)
	{
		Uint4 size = strlen(str);
		if (!current_block->check_space_size(size))
		{
			advance_block();
		}
		current_block->append_data(str, size);
		return *this;
	}
	ResultsBuffer& operator << (const string& str)
	{
		this->operator<<(str.c_str());
		return *this;
	}
	ResultsBuffer& operator << (Int8 n)
	{
		char buffer[64];
		Int8ToString(buffer, n);
		this->operator<<(buffer);
		return *this;
	}
	ResultsBuffer& operator << (Uint8 n)
	{
		char str[64];
		Uint8ToString(str, n);
		this->operator<<(str);
		return *this;
	}
	ResultsBuffer& operator << (Int4 n)
	{
		char buffer[64];
		Int8ToString(buffer, n);
		this->operator<<(buffer);
		return *this;
	}
	ResultsBuffer& operator << (Uint4 n)
	{
		char str[64];
		Uint8ToString(str, n);
		this->operator<<(str);
		return *this;
	}
	ResultsBuffer& Append(const char* str, Int4 size)
	{
		if (!current_block->check_space_size(size))
		{
			advance_block();
		}
		current_block->append_data(str, size);
		return *this;
	}
	void write_results(FILE* file)
	{
		flush(file);
	}
	
public:
	static void Uint8ToString(char str[], Uint8 n)
	{
		int num_digits = 0, index;
		Uint8 k = n;
		do
		{
			++num_digits;
			k /= 10;
		} while (k > 0);
		for (index = num_digits - 1; index >= 0; --index)
		{
			str[index] = (n % 10) + 48;
			n /= 10;
		}
		str[num_digits] = '\0';
	}

	static void Int8ToString(char str[], Int8 n)
	{
		int index = 0;
		if (n < 0)
		{
			str[index++] = '-';
			n = -n;
		}
		Uint8ToString(str + index, n);
	}
	
private:
	ResultBlk* blocks;
	ResultBlk* current_block;
};

typedef ResultsBuffer OutputBuffer;


#endif // RESULTS_BUFFER_H_

