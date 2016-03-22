/** @file index.h
 *  Implemente the FMD-INDEX used in bwa (https://github.com/lh3/bwa.git)
 *  We implement it here in C++ for convenience.
 */

#ifndef INDEX_H
#define	INDEX_H

#include "dbinfo.h"
#include "query_info.h"
#include "memallocator.h"

/* seed information */
struct MEM
{
	Int4 qoff;
	Int4 context;
	Int8 soff;
	Int8 diag;
	
	Int8 sid;
	Int4 block_id;
	Int4 block_offset;
};

/**
 * The bi-interval
 */
struct BwtIntv
{
    Uint8 x[3];
};

/**
 * The bi-interval is used during seeding phase.
 * To store the search results, only x[0] and x[2] in BwtIntv are needed.
 * In SearchInterval, k = x[0], l = x[2].
 */
struct SearchInterval
{
	Uint8 k,l;
	Int4 query_id, q_off;
};

/**
 * The FMD-INDEX, originally implemented in BWA.
 * nadb is the concatenation of the subjects in database, represented in blastna format.
 * ftable is the lookup table for databse
 */
class FMIndex
{
private:
    DbInfo* dbinfo;
    const char* dbname;
    char* src;
    
public:
    Uint8* sa;
    Uint4* bwt;
    Uint8  bwt_size;
    Uint8  seq_len;
    Uint8  L2[5];
    Uint8  primary;
    Uint1* pac;
    Uint1* nadb;
    Int4   lut_size;
    BwtIntv* ftable;
    Uint4  cnt_table[256];
    
public:
    FMIndex(const char* name, Uint8 r);
	FMIndex(const char* name);
    ~FMIndex();
    DbInfo* GetDbInfo() { return dbinfo; }
    void Print();
    
	// The following functions are used to build fm-index
	// To build fm-index, call BuildIndex()
    void BuildIndex();
    void RestoreBwt(const char* fn);
    void RestoreBwt2();
    void RestoreSa();
    void DumpBwt(const char* fn);
    void DumpBwt2(const char* fn);
    void UpdateBwt();
    void CalSA();
    void ConstructFTable();
    void MapBwt();
    void Destroy();
    void ConstructBwt();
	
private:
	void GenerateSuffixArrayFileName(char name[]);
	void GenerateBwtFileName(char name[]);
    
public:
    static const Uint8 kOccIntvShift = 7;
    static const Uint8 kOccInterval  = (1ull<<kOccIntvShift);
    static const Uint8 kOccIntvMask  = kOccInterval - 1;   
    static const Uint8 kSaIntv = 8;
    static const Int4  kLutSize = 12;
	static const char* index_build_name;
    
public:
	// The following functions are auxilliary tools for seeding
    void   Bwt2Occ(Uint8 k, Uint8 l, Uint1 base, Uint8& ok, Uint8& ol);
    int    BwtOccAux(Uint8 y, int c);
    Uint8  BwtOcc(Uint1 base, Uint8 pos);
    Uint4* GetBwtIntvAddr(Uint8 pos);
    Uint8  InvPsi(Uint8 pos);
    void   BwtOcc4(Uint8 k, Uint8 cnt[4]);
    void   Bwt2Occ4(Uint8 k, Uint8 l, Uint8 cntk[4], Uint8 cntl[4]);
    void   Extend(const BwtIntv* ik, BwtIntv ok[4], int is_back);
    Uint4 BackwardSearch(Uint1* seq, int len, Uint8& rk, Uint8& rl);
    Uint4 ForwardSearch(Uint1* seq, int len, Uint8& rk, Uint8& rl);
    
public:
	// Find the bi-intervals of seeds.
    Int4 Seeding(Uint1* query, Int4 query_length, Int4 word_size,
                 Int4 q_start, Int4 context, 
			     cy_utility::SimpleArray<SearchInterval>& intvs);
	// Return sa[k]
    Uint8 BwtSa(Uint8 k);
	// For each bi-interval [x1,x2,x3], and each k in [x1,x2,x3]
	// find out sa[k]
    Int4 LocateSeeds(cy_utility::SimpleArray<SearchInterval>& sis, 
					 cy_utility::SimpleArray<MEM*>& seeds, 
				     Int8 dblen, Int4 word_size, QueryInfo* query_info,
					 SmallObjAllocator& soa) ;
	Int4 LocateSeeds(cy_utility::SimpleArray<SearchInterval>& sis, 
					 cy_utility::SimpleArray<MEM>& seeds, 
				     Int8 dblen, Int4 word_size, QueryInfo* query_info) ;
	Uint8 BwtSa1(Uint8 k) { return sa[k]; }
};

// requirement: (OCC_INTERVAL % 16 == 0)
#define OCC_INTERVAL 0x80

#define bwt_bwt(k) (bwt[((k)>>7<<4) + sizeof(Uint8) + (((k)&0x7f)>>4)])

#define bwt_B0(k) (bwt_bwt((k))>>((~(k)&0xf)<<1)&3)
#define bwt_occ_intv(k) (bwt + ((k)>>7<<4))

#define bwt_set_intv(c, ik) ((ik).x[0] = L2[(int)(c)]+1, (ik).x[2] = L2[(int)(c)+1]-L2[(int)(c)], (ik).x[1] = L2[3-(c)]+1)

inline int
FMIndex::BwtOccAux(Uint8 y, int c)
{
    y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
    y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
    return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

inline Uint8
FMIndex::BwtOcc(Uint1 base, Uint8 pos)
{
    Uint8 n;
    Uint4 *p, *end;

    if (pos == (Uint8)seq_len) return L2[base + 1] - L2[base];
    if (pos == (Uint8)(-1)) return 0;
    if (pos >= primary) --pos;

    n = ((Uint8*)(p = GetBwtIntvAddr(pos)))[base];
    p += sizeof(Uint8);

    end = p + (((pos>>5) - ((pos&~kOccIntvMask)>>5))<<1);
    for (; p < end; p += 2) n += BwtOccAux((Uint8)p[0]<<32 | p[1], base);

    n += BwtOccAux(((Uint8)p[0]<<32 | p[1]) &~ ((1ull<<((~pos&31)<<1)) - 1), base);
    if (base == 0) n -= ~pos&31;

    return n;
}

inline Uint4*
FMIndex::GetBwtIntvAddr(Uint8 pos)
{
    Uint8 offset = (pos>>7)<<4;
    return bwt + offset;
}

inline Uint8 FMIndex::BwtSa(Uint8 k)
{
	Uint8 s = 0, mask = kSaIntv - 1;
	while (k & mask)
	{
		++s;
		k = InvPsi(k);
	}

	return s + sa[k / kSaIntv];
}

inline Uint8 FMIndex::InvPsi(Uint8 pos)
{
	Uint8 x = pos - (pos > primary);
	x = bwt_B0(x);
	x = L2[x] + BwtOcc(x, pos);
	return pos == primary ? 0 : x;
}

#endif	/* INDEX_H */

