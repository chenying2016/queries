#ifndef DEF_H
#define	DEF_H

#include <stdint.h>

typedef int8_t    Int1;
typedef uint8_t   Uint1;
typedef int16_t   Int2;
typedef uint16_t  Uint2;
typedef int32_t   Int4;
typedef uint32_t  Uint4;
typedef int64_t   Int8;
typedef uint64_t  Uint8;

typedef Int8 index_t;

#ifdef bool
typedef bool      Boolean;
#define TRUE      true
#define FALSE     false
#else
typedef Int1      Boolean;
#define TRUE      1
#define FALSE     0
#endif

#ifndef NULL
#define NULL 0
#endif

#include <cassert>
#define ASSERT    assert

#define ABS(c)    ((c)>=0 ? (c) : (-(c)))
#define ABS_BIT_CLEAR(c) (((c)<<1)>>1)

#ifndef MIN
#define MIN(a,b)  ((a)>(b) ? (b) : (a))
#endif

#ifndef MAX
#define MAX(a,b)  ((a)>(b) ? (a) : (b))
#endif

#ifndef INT4_MAX
#define INT4_MAX   2147483647
#endif

#ifndef UINT4_MAX
#define UINT4_MAX  4294967295U
#endif

#ifndef INT4_MIN
#define INT4_MIN  (-2147483647-1)
#endif

#ifndef INT2_MAX
#define INT2_MAX  32767
#endif

#ifndef INT2_MIN
#define INT2_MIN (-32768)
#endif

#ifndef INT1_MAX
#define INT1_MAX 127
#endif

#ifndef INT1_MIN
#define INT1_MIN (-128)
#endif

#ifndef UINT1_MAX
#define UINT1_MAX 256
#endif

#ifndef DIM
#define DIM(A)     (sizeof(A)/sizeof((A)[0]))
#endif

#ifndef INT4_MAX
#define INT4_MAX   2147483647
#endif

#ifndef UINT4_MAX
#define UINT4_MAX  4294967295U
#endif

#ifndef INT4_MIN
#define INT4_MIN  (-2147483647-1)
#endif

#ifndef INT2_MAX
#define INT2_MAX  32767
#endif

#ifndef INT2_MIN
#define INT2_MIN (-32768)
#endif

#ifndef QUERY_DELIMITER
#define QUERY_DELIMITER 0xf
#define CHAR_QUERY_DELIMITER '-'
#endif

#ifndef _set_pac
#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l)    ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
#endif

#define EXTRACTNA2BASE(pac, l) _get_pac(pac, (l))
#define SETNA2BASE(pac, l, c)  _set_pac(pac, (l), (c))

#define BLAST2NA_SIZE 4     /**< Size of compressed nucleic acid alphabet */
#define BLASTNA_SIZE 16     /**< Size of nucleic acid alphabet */
#define BLASTAA_SIZE 28     /**< Size of aminoacid alphabet */


#define BLASTNA_SEQ_CODE 99 /**< Identifies the blastna alphabet, for use in  blast only. */
#define BLASTAA_SEQ_CODE 11 /**< == Seq_code_ncbistdaa */
#define NCBI4NA_SEQ_CODE 4  /**< == Seq_code_ncbi4na */

/** Natural log(2) */
#define NCBIMATH_LN2    0.69314718055994530941723212145818

#define MASK64BIT_L32 (((Uint8(1)) << 32) - 1)
#define MASK64BIT_H32 (MASK64BIT_L32 << 32);

template <typename T, Int4 dim>
struct Point
{
private:
    T data[dim];
public:
    T& operator[](Int4 index) { return data[index]; }
    const T& operator[] (Int4 index) const { return data[index]; }

    Point() {}
    Point(const Point& pt)
    {
        for (Int4 i = 0; i < dim; ++i) data[i] = pt[i];
    }
    Point& operator=(const Point& pt)
    {
        for (Int4 i = 0; i < dim; ++i) data[i] = pt[i];
        return *this;
    }
    ~Point() {}
};

typedef Point<Int4, 2>    Point2d_Int4;
typedef Point<Uint4, 2>   Point2d_Uint4;
typedef Point<Int8, 2>    Point2d_Int8;
typedef Point<Uint8, 2>   Point2d_Uint8;

#define __sfree(x) { free(*x); *x = NULL; }
#define sfree(x) { void** __sfree_t = (void**)(void*)&(x); __sfree(__sfree_t); }

/**
 * A macro expression that returns 1, 0, -1 if a is greater than,
 * equal to or less than b, respectively.  This macro evaluates its
 * arguments more than once.
 */
#ifndef BLAST_CMP
#define BLAST_CMP(a,b) ((a)>(b) ? 1 : ((a)<(b) ? -1 : 0))
#endif

inline
Uint1 Blastna2Na2(Uint1 c)
{
    const  Uint1 code_table[16] = {0,1,2,3,0,0,0,0,0,0,0,0,0,0,0,QUERY_DELIMITER};
    return code_table[c];
}

inline Int4 Gcd(Int4 a, Int4 b)
{
	Int4   c;

	b = ABS(b);
	if (b > a)
		c=a, a=b, b=c;

	while (b != 0) {
		c = a%b;
		a = b;
		b = c;
	}
	return a;
}

inline Int4 
Gdb3(Int4* a, Int4* b, Int4* c)
{
    Int4 g;
    if (*b == 0) 
        g = Gcd(*a, *c);
    else 
        g = Gcd(*a, Gcd(*b, *c));
    if (g > 1) {
        *a /= g;
        *b /= g;
        *c /= g;
    }
    return g;
}

#define print_msg(out, msg) \
    do { \
        (out) << "[" << __FILE__ << ", " << __func__ << ", " << __LINE__ << "] " \
              << (msg) << "\n"; \
    } while (0) 

#define error_and_exit(msg) \
    do { \
        print_msg(std::cerr, msg); \
        exit(1); \
    } while (0)

#define warning(msg) \
    do { \
        std::cerr << "[" << __FILE__ << ", " << __func__ << ", " << __LINE__ << "] warning: " \
                  << msg << "\n"; \
    } while (0)

#endif	/* DEF_H */

