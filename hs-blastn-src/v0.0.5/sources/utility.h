#ifndef ULTILITY_H_
#define ULTILITY_H_ 1

#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <stdint.h>
#include <stdarg.h>
#include <sys/time.h>

// name space cy_utility
// Implement some simple datastructures
// Log, for tracing, print error massages
// NString, integer <=> string
// MemoryAllocator, a wrapper for memory allocator, including some error handling
// FileOperator, open, close, read and write files
// SimpleArray, an light-weight, std::vector-like dynamic array.
// Timer, records runtime.

namespace cy_utility
{
	extern const char* kAligner;
	
    struct Types
    {
        typedef int8_t   Int1;
        typedef uint8_t  Uint1;
        typedef int32_t  Int4;
        typedef uint32_t Uint4;
        typedef int64_t  Int8;
        typedef uint64_t Uint8;

        typedef Uint8    SIZE_TYPE;
    };

    struct Log
    {
#define ErrorMassageToString \
            const int kBufferSize = 2048; \
            char buf[kBufferSize]; \
            int msg_size; \
            va_list ap; \
            va_start(ap, format); \
            msg_size = vsprintf(buf, format, ap); \
            va_end(ap); \
            assert(msg_size < kBufferSize); \
            using std::clog; \
            using std::cerr; \
            using std::endl;
	
		static const int kDoTrace = 0;
		static int DoTrace() 
		{
			return kDoTrace;
		}
		static void Trace(const char* class_name, const char* format, ...)
		{
			if (!DoTrace()) return;
			
			ErrorMassageToString;
			if (class_name)
				clog << "[" << class_name << "] ";
			clog << buf << endl << endl;
		}

        static void ErrorAndExit(const char* class_name, const char* format, ...)
        {
            ErrorMassageToString;

            if (class_name != NULL)
                cerr << "[" << class_name << "] ";
            cerr << "Error: "
                 << buf << endl;
            exit(1);
        }

        static void Warning(const char* class_name, const char* format, ...)
        {
            ErrorMassageToString;

            if (class_name != NULL)
                cerr << "[" << class_name << "] ";
            cerr << "Warning: "
                 << buf << endl;
        }

        static void LogMsg(const char* class_name, const char* format, ...)
        {
            ErrorMassageToString;

            if (class_name != NULL)
                clog << "[" << class_name << "] ";
            clog << buf << endl;
        }
    };

	/// Number to string conversion flags.
    ///
    /// NOTE: 
    ///   If specified base in the *ToString() methods is not default 10,
    ///   that some flags like fWithSign and fWithCommas will be ignored.
    enum ENumToStringFlags {
        fWithSign        = (1 <<  6), ///< Prefix the output value with a sign
        fWithCommas      = (1 <<  7), ///< Use commas as thousands separator
        fDoubleFixed     = (1 <<  8), ///< Use n.nnnn format for double
        fDoubleScientific= (1 <<  9), ///< Use scientific format for double
        fDoublePosix     = (1 << 10), ///< Use C locale
        fDoubleGeneral   = fDoubleFixed | fDoubleScientific,
        fDS_Binary               = (1 << 11),
        fDS_NoDecimalPoint       = (1 << 12),
        fDS_PutSpaceBeforeSuffix = (1 << 13),
        fDS_ShortSuffix          = (1 << 14),
        fDS_PutBSuffixToo        = (1 << 15)
    };
    typedef int TNumToStringFlags;    ///< Bitwise OR of "ENumToStringFlags"
	// A maximal double precision used in the double to string conversion
#if defined(NCBI_OS_MSWIN)
    const int kMaxDoublePrecision = 200;
#else
    const int kMaxDoublePrecision = 308;
#endif
// A maximal size of a double value in a string form.
// Exponent size + sign + dot + ending '\0' + max.precision
const int kMaxDoubleStringSize = 308 + 3 + kMaxDoublePrecision;
	
    struct NString
    {
        typedef Types::SIZE_TYPE SIZE_TYPE;
		

        static int CheckInteger(const char* str)
        {
            while (*str != '\0' && !isspace(*str))
            {
                if (!isdigit(*str)) return 0;
                ++str;
            }

            return 1;
        }

        static int GetBitsOfInteger(SIZE_TYPE n)
        {
            int r = 0;
            do
            {
                ++r;
                n /= 10;
            } while (n > 0);
            return r;
        }

        static int StringToInteger(const char* str, SIZE_TYPE& n)
        {
            if (!CheckInteger(str))
                return 1;
            n = 0;
            while (*str != '\0' && !isspace(*str))
            {
                n = n * 10 + *str - 48;
                ++str;
            }

            return 0;
        }

        static void IntegerToString(char str[], SIZE_TYPE n)
        {
            int bits = GetBitsOfInteger(n);
			str[bits] = '\0';
            int i;
            for (i = 0; i < bits; ++i)
            {
                str[bits - i - 1] = (n % 10) + 48;
                n /= 10;
            }
        }
		
		static int DoubleToString(double value, int precision, 
								  char* buf, int buf_size,
								  TNumToStringFlags flags)
		{
			char buffer[kMaxDoubleStringSize];
			int n = 0;
			if ((flags & fDoublePosix) && (!finite(value) || value == 0.))
			{
				if (value == 0.)
				{
					double zero = 0.;
					if (memcpy(&value, &zero, sizeof(double)) == 0)
					{
						strcpy(buffer, "0");
						n = 2;
					}
					else
					{
						strcpy(buffer, "-0");
						n = 3;
					}
				}
				else if (isnan(value))
				{
					strcpy(buffer, "NaN");
					n = 4;
				}
				else if (value > 0.)
				{
					strcpy(buffer, "INF");
					n = 4;
				}
				else
				{
					strcpy(buffer, "-INF");
					n = 5;
				}
			}
			else
			{
				if (precision > kMaxDoublePrecision)
					precision = kMaxDoublePrecision;
				const char* format;
				switch (flags & fDoubleGeneral)
				{
					case fDoubleScientific:
						format = "%.*e";
						break;
					case fDoubleGeneral:
						format = "%.*g";
						break;
					case fDoubleFixed:
					default:
						format = "%.*f";
						break;
				}
				n = sprintf(buffer, format, precision, value);
				if (n < 0)
					n = 0;
				if (flags & fDoublePosix)
				{
					struct lconv* conv = localeconv();
					if ('.' != *(conv->decimal_point))
					{
						char* pos = strchr(buffer, *(conv->decimal_point));
						if (pos) *pos = '.';
					}
				}
			}
			int n_copy = std::min(n, buf_size);
			memcpy(buf, buffer, n_copy);
			return n_copy;
		}
    };
	
    struct MemoryAllocator
    {
        typedef Types::SIZE_TYPE SIZE_TYPE;
        static const char* kClassName;

        static void* __malloc(SIZE_TYPE size)
        {
            void* p = malloc(size);
            if (p == NULL)
                Log::ErrorAndExit(kClassName, "memory allocation failed.");
            return p;
        }
		
		static void* __calloc(SIZE_TYPE size)
		{
			void* p = calloc(1, size);
			if (p == NULL)
				Log::ErrorAndExit(kClassName, "memory allocation failed.");
			return p;
		}

        static void* __realloc(void* p, SIZE_TYPE new_size)
        {
            void* r = realloc(p, new_size);
            if (r == NULL)
                Log::ErrorAndExit(kClassName, "memory reallocation failed.");
            return r;
        }

        static void* __free(void* p)
        {
            if (p != NULL) free(p);
            return NULL;
        }
    };

    struct FileOperator
    {
        typedef Types::SIZE_TYPE SIZE_TYPE;
        static const SIZE_TYPE kSegmentSize = 200 * (1ULL << 20);

        static FILE* openfile(const char* class_name, const char* file_name, const char* flag)
        {
            if (file_name == NULL)
                Log::ErrorAndExit(class_name, "file name is empty.");
            FILE* file = fopen(file_name, flag);
            if (file == NULL)
                Log::ErrorAndExit(class_name, "cannot open file %s for reading.", file_name);
            return file;
        }

        static FILE* closefile(FILE* file)
        {
            if (file == NULL) return NULL;
            fclose(file);
            return NULL;
        }

        static void read_file(const char* class_name, const char* file_name, FILE* file, void* buf, SIZE_TYPE size)
        {
           SIZE_TYPE segment_size, left = size, offset = 0, read_size;
           char* p = (char*)buf;
           while (left > 0)
           {
               segment_size = kSegmentSize;
               if (segment_size > left) segment_size = left;
               read_size = fread(p + offset, 1, segment_size, file);
               if (read_size != segment_size)
               {
                   char err_msg[2048];
                   strcpy(err_msg, "while reading file ");
                   strcat(err_msg, file_name);
                   strcat(err_msg, "unexpected end of file.");

                   Log::ErrorAndExit(class_name, "%s", err_msg);
               }
               left -= segment_size;
               offset += segment_size;
           }
        }

        static void write_file(const char* class_name, const char* file_name, FILE* file, const void* buf, SIZE_TYPE size)
        {
            SIZE_TYPE segment_size, left = size, offset = 0, write_size;
            char* p = (char*)buf;
            while (left > 0)
            {
                segment_size = kSegmentSize;
                if (segment_size > left) segment_size = left;
                write_size = fwrite(p + offset, 1, segment_size, file);
                if (write_size != segment_size)
                {
                    char err_msg[2048];
                    strcpy(err_msg, "error occur when writing file ");
                    strcat(err_msg, file_name);
                    Log::ErrorAndExit(class_name, "%s", err_msg);
                }
                left -= segment_size;
                offset += segment_size;
            }
        }
    };

    template <typename T>
    struct SimpleArray
    {
        typedef Types::SIZE_TYPE SIZE_TYPE;
		static const SIZE_TYPE kMinArrSize = 64;

        T* arr;
        SIZE_TYPE arr_size;
        SIZE_TYPE alloc_size;

        SimpleArray(SIZE_TYPE n = kMinArrSize)
        {
			if (n < kMinArrSize) n = kMinArrSize;
            alloc_size = n;
            arr_size = 0;
            arr = (T*)MemoryAllocator::__malloc(sizeof(T) * alloc_size);
        }
        void destroy()
        {
            arr = (T*)MemoryAllocator::__free(arr);
            arr_size = alloc_size = 0;
        }
        ~SimpleArray()
        {
            destroy();
        }
        void reserve(SIZE_TYPE new_size)
        {
            if (new_size <= alloc_size) return;

            arr = (T*)MemoryAllocator::__realloc(arr, sizeof(T) * new_size);
            alloc_size = new_size;
        }
		void set_array_size(SIZE_TYPE size)
		{
			assert(size <= arr_size);
			arr_size = size;
		}
        void* get_data()
        {
            return arr;
        }
		void set_data(T* src, SIZE_TYPE array_size, SIZE_TYPE mem_size, int no_free = true)
		{
			if (!no_free)
				destroy();
			arr = src;
			arr_size = array_size;
			alloc_size = mem_size;
		}
        void push_back(const T* src, SIZE_TYPE size)
        {
            if (arr_size + size > alloc_size)
            {
                while (arr_size + size > alloc_size) 
                    alloc_size *= 2;

                arr = (T*)MemoryAllocator::__realloc(arr, sizeof(T) * alloc_size);
            }

            memcpy(arr + arr_size, src, sizeof(T) * size);
            arr_size += size;
        }
        void push_back(T item)
        {
            if (arr_size >= alloc_size)
            {
                alloc_size *= 2;
                arr = (T*)MemoryAllocator::__realloc(arr, sizeof(T) * alloc_size);
            }
            arr[arr_size++] = item;
        }
        void push_back_simple(T item)
        {
            arr[arr_size++] = item;
        }
        T& operator[](SIZE_TYPE index)
        {
            return arr[index];
        }
        const T& operator[](SIZE_TYPE index) const
        {
            return arr[index];
        }
        SIZE_TYPE size()
        {
            return arr_size;
        }
        SIZE_TYPE size() const
        {
            return arr_size;
        }
        void clear()
        {
            arr_size = 0;
        }
        void clear() const
        {
            arr_size = 0;
        }
    };

    struct Timer
    {
        struct timeval m_Start, m_End;
        void start()
        {
            gettimeofday(&m_Start, NULL);
        }
        void end()
        {
            gettimeofday(&m_End, NULL);
        }
        double get_elapsed_time()
        {
            double r;
            r = 1.0 * (m_End.tv_sec - m_Start.tv_sec) + 1.0 * (m_End.tv_usec - m_Start.tv_usec) / 1000000;
            return r;
        }
    };
}

#endif // ULTILITY_H_
