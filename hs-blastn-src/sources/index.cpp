#include "index.h"
#include <vector>
#include <iostream>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

const char* FMIndex::index_build_name = "IndexBuilder";
static const char* kSuffixArrayFileSuffix = ".sa";
static const char* kBwtFileSuffix = ".bwt";

void FMIndex::GenerateSuffixArrayFileName(char name[])
{
	if (dbname == NULL)
	{
		cy_utility::Log::ErrorAndExit(index_build_name, "database name is empty.");
	}
	strcpy(name, dbname);
	strcat(name, kSuffixArrayFileSuffix);
}

void FMIndex::GenerateBwtFileName(char name[])
{
	if (dbname == NULL)
	{
		cy_utility::Log::ErrorAndExit(index_build_name, "database name is empty.");
	}
	strcpy(name, dbname);
	strcat(name, kBwtFileSuffix);
}

void FMIndex::Print()
{	
	using std::cerr;
	using std::endl;
	
	cerr << index_build_name << endl;
	cerr << "\tprimary = " << primary << endl;
	cerr << "\tnumber of letters = " << seq_len << endl;
	cerr << "\tbwt_size = " << bwt_size << endl;
	cerr << "\tlut_size = " << lut_size << endl;
	
	cerr << "\t";
	for (int i = 0; i < 5; ++i)
		cerr << "L[" << i << "] = " << L2[i] << "\t";
	cerr << endl;
}

void FMIndex::CalSA()
{
	using namespace cy_utility;
	Timer timer;
	std::cerr << std::endl;
	Log::LogMsg(index_build_name, "construct suffix array.");
	timer.start();
	
	using std::cerr;
	using std::endl;
	
	Uint8 n_sa = (seq_len + kSaIntv) / kSaIntv;
	sa = (Uint8*)cy_utility::MemoryAllocator::__malloc(n_sa * sizeof(Uint8));
	
	Uint8 s = seq_len, rs = 0;
	Uint8 i;
	Uint8 report_size = 100000000, size, left = seq_len + 1, processed = 0;
	
	while (left > 0)
	{
		size = std::min(report_size, left);
		for (i = 0; i < size; ++i)
		{
			if (rs % kSaIntv == 0) sa[rs / kSaIntv] = s;
			--s;
			rs = InvPsi(rs);
		}
		left -= size;
		processed += size;
		cerr << "[" << index_build_name << "] "
			 << processed << "/" << seq_len + 1 << " suffix array elements processed.\n";
	}
	
	if (rs % kSaIntv == 0) sa[rs / kSaIntv] = s;
	sa[0] = (Uint8)(-1);
    
    char name[1024];
	GenerateSuffixArrayFileName(name);
	FILE* file = cy_utility::FileOperator::openfile(index_build_name, name, "w");
	cy_utility::FileOperator::write_file(index_build_name, name, file, sa, n_sa * sizeof(Uint8));
    sa = (Uint8*)cy_utility::MemoryAllocator::__free(sa);
	file = cy_utility::FileOperator::closefile(file);
	
	timer.end();
	double dur = timer.get_elapsed_time();
	Log::LogMsg(index_build_name, "done. Time elapsed: %f secs.", dur);
}

Int4 FMIndex::LocateSeeds(cy_utility::SimpleArray<SearchInterval>& sis, 
						  cy_utility::SimpleArray<MEM>& seeds, 
					      Int8 dblen, Int4 word_size, 
						  QueryInfo* query_info) {
	SearchInterval* intvs = (SearchInterval*) sis.get_data();
	Int4 num_sis = sis.size();
	if (num_sis == 0) return 0;

	Int4 num_hits = 0;
	MEM m; m.sid = 0; m.diag = 0; m.block_id = m.block_offset = 0; // g++ complains if these quantities are not initialized
	Int8 soff;
	Uint8 k;
	Int4 i;
	
	cy_utility::Log::Trace(index_build_name, "Number of intervals: %d", num_sis);
	
	Int8 t = dblen;

	for (i = 0; i < num_sis; ++i) {
		SearchInterval& intv = intvs[i];
		
		for (k = 0; k < intv.l; ++k) {             
			soff = BwtSa(k + intv.k);
			
			if (soff >= t)
			{
				m.context = intv.query_id * 2 + 1;
				m.soff = 2 * dblen - word_size - soff;
				m.qoff = query_info->GetSeqLength(m.context) - word_size - intv.q_off;
			}
			else
			{
				m.context = intv.query_id * 2;
				m.soff = soff;
				m.qoff = intv.q_off;
			}   
			
			seeds.push_back(m);
		}
	}

	return num_hits;
}

void FMIndex::ConstructBwt()
{
    ASSERT(sa != NULL);
    ASSERT(sa[0] == seq_len);
    
    bwt_size = (seq_len + 15) >> 4;
    bwt = (Uint4*)calloc(sizeof(Uint4), bwt_size);
    ASSERT(bwt != NULL);
    
    Uint8 i, j = 0;
    Uint1* buf = (Uint1*)calloc(sizeof(Uint1), seq_len);
    
    for (i = 0; i <= seq_len; ++i)
    {
        if (sa[i] == 0) primary = i;
        else buf[j++] = pac[sa[i] - 1];
    }
    ASSERT(j == seq_len);
    
    for (i = 0; i < seq_len; ++i)
        bwt[i >> 4] |= buf[i] << ((15 - (i & 15)) << 1);
    
    free(buf);
}

FMIndex::FMIndex(const char* name)
{
    dbname = name;
    dbinfo = new DbInfo(name);
	pac = NULL;
	nadb = NULL;
    src = NULL;
    bwt = NULL;
    ftable = NULL;
    bwt_size = 0;
    seq_len = 0;
    primary = 0;
    lut_size = kLutSize;
    memset(L2, 0, sizeof(Uint8) * 5);
}

FMIndex::~FMIndex()
{
    Destroy();
    
    if (sa) free(sa); sa = NULL;
    
    if (dbinfo != NULL) delete dbinfo;
}

void FMIndex::RestoreSa()
{
	using std::cerr;
	using std::endl;
    using std::clog;
	
	char name[2048];
	GenerateSuffixArrayFileName(name);
	FILE* file = cy_utility::FileOperator::openfile(index_build_name, name, "r");
	fseek(file, 0ULL, SEEK_END);
	Uint8 size = ftell(file);
	fseek(file, 0ULL, SEEK_SET);

    Uint8 GB = 1;
    GB = GB << 30;
    double gb = 1.0 * size / GB;
    //clog << "\tLoading " << name << ", size = " << gb << " GB\n";
    fprintf(stderr, "\tLoading %s, size = %.1gGB\n", name, gb);

	sa = (Uint8*)cy_utility::MemoryAllocator::__malloc(size);
	cy_utility::FileOperator::read_file(index_build_name, name, file, sa, size);
	file = cy_utility::FileOperator::closefile(file);
}

void FMIndex::Destroy()
{   
    bwt = NULL;
    ftable = NULL;
    bwt_size = 0;
    seq_len = 0;
    primary = 0;
    memset(L2, 0, sizeof(Uint8) * 5);
	
	if (src != NULL)
	{
		free(src);
		src = NULL;
	} else 
	{
		if (bwt)
		{
			free(bwt); bwt = NULL;
		}
		if (ftable)
		{
			free(ftable); ftable = NULL;
		}
		if (sa)
		{
			free(sa); sa = NULL;
		}
		if (pac)
		{
			free(pac); pac = NULL;
		}
		if (nadb)
		{
			free(nadb); nadb = NULL;
		}
	}
}

void bwt_gen_cnt_table(Uint4 cnt_table[256])
{
    int i, j;
    for (i = 0; i != 256; ++i)
    {
        Uint4 x = 0;
        for (j = 0; j != 4; ++j)
            x |= (((i & 3) == j) + ((i >> 2 & 3) == j) + ((i >> 4 & 3) == j) + (i >> 6 == j)) << (j << 3);
        cnt_table[i] = x;
    }
}

#include "bwt_gen.h"

void FMIndex::BuildIndex()
{   

    dbinfo->MakePostedDate();
    dbinfo->BuildDbInfo();

	using namespace cy_utility;
	Timer timer;
	double dur;
	std::clog << std::endl;
	Log::LogMsg(index_build_name, "packing database.");
	timer.start();
	
    const Uint1 table[16] = {0,1,2,3,0,0,0,0,0,0,0,0,0,0,0,0};
    Uint8 dblen = dbinfo->GetDbLength();
    Uint8 pac_bytes = (2 * dblen + 3) / 4;
    pac = (Uint1*)calloc(sizeof(Uint1), pac_bytes);
    ASSERT(pac != NULL);
    Uint1* nadb = (Uint1*)dbinfo->GetDb();
    Uint8 i, j;
    Uint1 c;
    seq_len = 2 * dblen;
    
    for (i = 0; i < dblen; ++i)
    {
        c = table[nadb[i]];
        _set_pac(pac, i, c);
    }
    
    for (i = 0; i < dblen; ++i)
    {
        c = 3 - table[nadb[dblen - i - 1]];
        j = dblen + i;
        _set_pac(pac, j, c);
    }
	
	timer.end();
	dur = timer.get_elapsed_time();
	Log::LogMsg(index_build_name, "done. Time eaplsed: %f secs.", dur);
    
    char name[1024];
    char name2[1024];
    strcpy(name, dbname);
    strcat(name, ".sa");
    strcpy(name2, dbname);
    strcat(name2, ".bwt");
    FILE* file = fopen(name, "w");
    ASSERT(file != NULL);
    i = fwrite(pac, 1, pac_bytes, file);
    ASSERT(i == pac_bytes);
    
    if (seq_len % 4 == 0)
    {
        c = 0;
        i = fwrite(&c, 1, 1, file);
        ASSERT(i == 1);
    }
    
    c = seq_len % 4;
    i = fwrite(&c, 1, 1, file);
    ASSERT(i == 1);
    fclose(file);
    
    dbinfo->Destroy();
    free(pac);
    pac = NULL;
    
	std::clog << std::endl;
	Log::LogMsg(index_build_name, "construct bwt.");
	timer.start();
    bwt_bwtgen(name, name2);
	timer.end();
	dur = timer.get_elapsed_time();
	Log::LogMsg(index_build_name, "done. Time elapsed: %f secs.", dur);
 
    RestoreBwt(name2);
    
    UpdateBwt();
    
    bwt_gen_cnt_table(cnt_table);
    
    ConstructFTable();
    
    DumpBwt2(name2);
    
	CalSA();

	free(bwt); bwt = NULL;
	free(ftable); ftable = NULL;
}

void FMIndex::RestoreBwt(const char* fn)
{	
	using namespace cy_utility;
	FILE* file = FileOperator::openfile(index_build_name, fn, "r");
	fseek(file, 0ULL, SEEK_END);
	bwt_size = (ftell(file) - sizeof(Uint8) * 5) >> 2;
	fseek(file, 0ULL, SEEK_SET);
	bwt = (Uint4*)MemoryAllocator::__malloc(bwt_size * sizeof(Uint4));
	
	FileOperator::read_file(index_build_name, fn, file, &primary, sizeof(Uint8));
	FileOperator::read_file(index_build_name, fn, file, L2 + 1, 4 * sizeof(Uint8));
	FileOperator::read_file(index_build_name, fn, file, bwt, bwt_size * sizeof(Uint4));
	seq_len = L2[4];
	
	file = FileOperator::closefile(file);
}

void FMIndex::RestoreBwt2()
{

    using std::clog;

	bwt_gen_cnt_table(cnt_table);
	
	if (dbname == NULL)
		cy_utility::Log::ErrorAndExit(index_build_name, "database file name is empty.");
	
	char name[2048];
	GenerateBwtFileName(name);
	FILE* bwt_file = cy_utility::FileOperator::openfile(index_build_name, name, "r");
	fseek(bwt_file, 0ULL, SEEK_END);
	Uint8 file_size = ftell(bwt_file);
	fseek(bwt_file, 0ULL, SEEK_SET);

    Uint8 GB = 1;
    GB = GB << 30;
    double gb = 1.0 * file_size / GB;
    fprintf(stderr, "\tLoading %s, size = %.1gGB\n", name, gb);

	src = (char*)cy_utility::MemoryAllocator::__malloc(file_size);
	cy_utility::FileOperator::read_file(index_build_name, name, bwt_file, src, file_size);
	bwt_file = cy_utility::FileOperator::closefile(bwt_file);
    
    Int8 index = 0;
    
    /// 1)
    primary = *(Uint8*)src;
    index += sizeof(Uint8);
    
    /// 2)
    seq_len = *(Uint8*)(src + index);
    index += sizeof(Uint8);
    
    /// 3)
    bwt_size = *(Uint8*)(src + index);
    index += sizeof(Uint8);
    
    /// 4)
    lut_size = *(Int4*)(src + index);
    index += sizeof(Int4);
    
    /// 5)
    memcpy(L2, src + index, sizeof(Uint8) * 5);
    index += sizeof(Uint8) * 5;
    
    /// 6)
    bwt = (Uint4*)(src + index);
    index += sizeof(Uint4) * bwt_size;
    
    /// 7)
    ftable = (BwtIntv*)(src + index);  
}

/*
 * Index Content
 * 1) primary  Uint8
 * 2) seq_len  Uint8
 * 3) bwt_size Uint8
 * 4) lut_size Int4
 * 5) L2       Uint8[5]
 * 6) bwt      Uint4*
 * 7) ftable   Uint4*
 */
void FMIndex::DumpBwt2(const char* fn)
{
	ASSERT(fn != NULL);
	FILE* file = cy_utility::FileOperator::openfile(index_build_name, fn, "w");
    
    Uint8 nws;

    /// 1)
	cy_utility::FileOperator::write_file(index_build_name, fn, file, &primary, sizeof(Uint8));
    
    /// 2)
	cy_utility::FileOperator::write_file(index_build_name, fn, file, &seq_len, sizeof(Uint8));
    
    /// 3)
	cy_utility::FileOperator::write_file(index_build_name, fn, file, &bwt_size, sizeof(Uint8));
    
    /// 4)
	cy_utility::FileOperator::write_file(index_build_name, fn, file, &lut_size, sizeof(Int4));
    
    /// 5)
	cy_utility::FileOperator::write_file(index_build_name, fn, file, L2, 5 * sizeof(Uint8));
    
    /// 6)
	cy_utility::FileOperator::write_file(index_build_name, fn, file, bwt, bwt_size * sizeof(Uint4));
    
    /// 7)
	Uint8 ftsize = (1ULL << (kLutSize << 1));
	cy_utility::FileOperator::write_file(index_build_name, fn, file, ftable, ftsize * sizeof(BwtIntv));
    
    file = cy_utility::FileOperator::closefile(file);
}

#define bwt_B00(bwt, k) ((bwt)[(k)>>4]>>((~(k)&0xf)<<1)&3)

void FMIndex::UpdateBwt()
{    
	std::cerr << std::endl;
	cy_utility::Timer timer;
	cy_utility::Log::LogMsg(index_build_name, "update bwt.");
	timer.start();
	
    Uint8 i, k, c[4], n_occ;
    Uint4* buf;
    
    n_occ = (seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
    bwt_size += n_occ * sizeof(Uint8);
	buf = (Uint4*)cy_utility::MemoryAllocator::__calloc(bwt_size * sizeof(Uint4));
    
    c[0] = c[1] = c[2] = c[3] = 0;
    for (i = k = 0; i < seq_len; ++i)
    {
        if (i % OCC_INTERVAL == 0)
        {
            memcpy(buf + k, c, sizeof(Uint8) * 4);
            k += sizeof(Uint8);
        }
        if (i % 16 == 0) buf[k++] = bwt[i/16];
        ++c[bwt_B00(bwt, i)];
    }
    memcpy(buf + k, c, sizeof(Uint8) * 4);
    ASSERT(k + sizeof(Uint8) == bwt_size);
    free(bwt);
    bwt = buf;
    
    // set L2
    L2[0] = 0;
    for (i = 1; i <=4; ++i)
    {
        L2[i] = L2[i-1] + c[i-1];
    }
	
	timer.end();
	double dur = timer.get_elapsed_time();
	cy_utility::Log::LogMsg(index_build_name, "done. Time elapsed: %f secs.", dur);
}

Uint4 FMIndex::BackwardSearch(Uint1* seq, int len, Uint8& rk, Uint8& rl)
{
    Uint8 k, l, ok, ol;
    int i;
    k = 0;
    l = seq_len;
    
    for (i = len - 1; i >= 0; --i)
    {
        Uint1 c = seq[i];
        Bwt2Occ(k - 1, l, c, ok, ol);
        k = L2[c] + ok + 1;
        l = L2[c] + ol;
        if (k > l) break;
    }
    if (k > l) return 0;
    rk = k;
    rl = l;
    return l - k + 1;
}

Uint4 FMIndex::ForwardSearch(Uint1* seq, int len, Uint8& rk, Uint8& rl)
{
    Uint8 k, l, ok, ol;
    int i;
    k = 0;
    l = seq_len;

    for (i = 0; i < len; ++i)
    {
        Uint1 c = seq[i];
        if (c > 3) return 0;
        Bwt2Occ(k - 1, l, c, ok, ol);
        k = L2[c] + ok + 1;
        l = L2[c] + ol;
        if (k > l) break;
    }
    if (k > l) return 0;
    rk = k;
    rl = l;
    return l - k + 1;
}

void FMIndex::Bwt2Occ(Uint8 k, Uint8 l, Uint1 base, Uint8& ok, Uint8& ol)
{
    uint64_t _k = k;
    uint64_t _l = l;
    if (k >= primary) --_k;
    if (l >= primary) --_l;

    if (_l / kOccInterval != _k / kOccInterval || k == (uint64_t) (-1) || l == (uint64_t) (-1))
    {
        ok = BwtOcc(base, k);
        ol = BwtOcc(base, l);
    }
    else
    {
        uint64_t m, n, i, j;
        uint32_t* p;
        if (k >= primary) --k;
        if (l >= primary) --l;

        n = ((uint64_t*) (p = GetBwtIntvAddr(k)))[base];
        p += sizeof (uint64_t);

        // compute ok
        j = (k >> 5) << 5;
        for (i = k / kOccInterval * kOccInterval; i < j; i += 32, p += 2)
            n += BwtOccAux(((uint64_t) p[0]) << 32 | p[1], base);
        m = n;
        n += BwtOccAux(((uint64_t) p[0] << 32 | p[1]) & ~((1ull << ((~k & 31) << 1)) - 1), base);
        if (base == 0) n -= ~k & 31;
        ok = n;

        // compute ol
        j = (l >> 5) << 5;
        for (; i < j; i += 32, p += 2)
            m += BwtOccAux(((uint64_t) p[0]) << 32 | p[1], base);
        m += BwtOccAux(((uint64_t) p[0] << 32 | p[1]) & ~((1ull << ((~l & 31) << 1)) - 1), base);
        if (base == 0) m -= ~l & 31;
        ol = m;
    }    
}

#define __occ_aux4(b)                                           \
	(cnt_table[(b)&0xff] + cnt_table[(b)>>8&0xff]		\
	 + cnt_table[(b)>>16&0xff] + cnt_table[(b)>>24])

void FMIndex::BwtOcc4(Uint8 k, Uint8 cnt[4])
{
    Uint8 x;
    Uint4 *p, tmp, *end;
    if (k == (Uint8)(-1))
    {
        memset(cnt, 0, sizeof(Uint8) * 4);
        return;
    }
    
    k -= (k >= primary);
    p = bwt_occ_intv(k);
    memcpy(cnt, p, sizeof(Uint8) * 4);
    p += sizeof(Uint8);
    end = p + ((k >> 4) - ((k&~kOccIntvMask) >> 4));
    for (x = 0; p < end; ++p) 
        x += __occ_aux4(*p);
    tmp = *p & ~((1U<<((~k&15)<<1))-1);
    x += __occ_aux4(tmp) - (~k&15);
    cnt[0] += x & 0xff;
    cnt[1] += x >> 8 & 0xff;
    cnt[2] += x >> 16 & 0xff;
    cnt[3] += x >> 24;
}

void FMIndex::Bwt2Occ4(Uint8 k, Uint8 l, Uint8 cntk[], Uint8 cntl[])
{
    Uint8 _k, _l;
    _k = k - (k >= primary);
    _l = l - (l >= primary);
    if (_l >> kOccIntvShift != _k >> kOccIntvShift || k == (Uint8) (-1) || l == (Uint8) (-1))
    {
        BwtOcc4(k, cntk);
        BwtOcc4(l, cntl);
    }
    else
    {
        Uint8 x, y;
        uint32_t *p, tmp, *endk, *endl;
        k -= (k >= primary); // because $ is not in bwt
        l -= (l >= primary);
        p = bwt_occ_intv(k);
        memcpy(cntk, p, 4 * sizeof (Uint8));
        p += sizeof (Uint8); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
        // prepare cntk[]
        endk = p + ((k >> 4) - ((k&~kOccIntvMask) >> 4));
        endl = p + ((l >> 4) - ((l&~kOccIntvMask) >> 4));
        for (x = 0; p < endk; ++p) x += __occ_aux4(*p);
        y = x;
        tmp = *p & ~((1U << ((~k & 15) << 1)) - 1);
        x += __occ_aux4(tmp) - (~k & 15);
        // calculate cntl[] and finalize cntk[]
        for (; p < endl; ++p) y += __occ_aux4(*p);
        tmp = *p & ~((1U << ((~l & 15) << 1)) - 1);
        y += __occ_aux4(tmp) - (~l & 15);
        memcpy(cntl, cntk, 4 * sizeof (Uint8));
        cntk[0] += x & 0xff;
        cntk[1] += x >> 8 & 0xff;
        cntk[2] += x >> 16 & 0xff;
        cntk[3] += x >> 24;
        cntl[0] += y & 0xff;
        cntl[1] += y >> 8 & 0xff;
        cntl[2] += y >> 16 & 0xff;
        cntl[3] += y >> 24;
    }    
}

void FMIndex::Extend(const BwtIntv* ik, BwtIntv ok[], int is_back)
{
    Uint8 tk[4], tl[4];
    int i;
    Bwt2Occ4(ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
    for (i = 0; i != 4; ++i)
    {
        ok[i].x[!is_back] = L2[i] + 1 + tk[i];
        ok[i].x[2] = tl[i] - tk[i];
    }
    ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= primary && ik->x[!is_back] + ik->x[2] - 1 >= primary);
    ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
    ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
    ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];    
}


void FMIndex::ConstructFTable()
{
#define __get_seed_index(c) \
	{ \
		offset = 0; \
		for(int k = 0; k < kLutSize - 1; ++k) \
		{ \
			assert(code[k] > 0); \
			offset = (offset<<2)|(code[k] - 1); \
		} \
		offset = (offset << 2) | c; \
	}
	
	cy_utility::Timer timer;
	std::cerr << std::endl;
	cy_utility::Log::LogMsg(index_build_name, "construct fast table");
	timer.start();
	
	Uint1 code[kLutSize];
	Uint1 kMaxCode = 3;
	Uint1 i;
	for (i = 0; i < kLutSize; ++i) code[i] = 0;
	
    Uint8 table_items = (1ULL << (kLutSize << 1));
    Uint8 table_bytes = sizeof(BwtIntv) * table_items;
	ftable = (BwtIntv*)cy_utility::MemoryAllocator::__malloc(table_bytes);
	BwtIntv ik, ok[4];
	BwtIntv warr[kLutSize];
	
	int level = 0;
	Uint4 offset;
	while (level >= 0)
	{
		if (code[level] > kMaxCode)
		{
			code[level] = 0;
			--level;
			continue;
		}
		
		if (level == 0)
		{
			bwt_set_intv(code[0], warr[0]);
			++code[0];
			++level;
			continue;
		}
		
		if (level == kLutSize - 1)
		{
			Extend(&warr[level - 1], ok, 0);
			for (i = 0; i <= kMaxCode; ++i)
			{
				__get_seed_index(i);
				ftable[offset] = ok[3 - i];
			}
			
			--level;
			continue;
		}
		
		Uint1 c = 3 - code[level];
		Extend(&warr[level - 1], ok, 0);
		warr[level] = ok[c];
		++code[level];
		++level;
	}
	
	timer.end();
	double dur = timer.get_elapsed_time();
	cy_utility::Log::LogMsg(index_build_name, "done. Time elapsed: %f secs.", dur);
}

static inline Uint4 GetSeedIndex(Uint1* s, Uint4 size)
{
    Uint4 i, ret = 0;
    for (i = 0; i < size; ++i)
    {
        ret = (ret<<2)|(*s++);
    }
    
    return ret;
}

struct EBwtIntv
{
    struct BwtIntv intv;
    Int4 ext_left;
};

#include <list>
#include <vector>

Int4 FMIndex::Seeding(Uint1* query, Int4 query_length, 
                      Int4 context, Int4 word_size, Int4 q_start,
                      cy_utility::SimpleArray<SearchInterval>& intvs)
{
    if (query_length < word_size) return 0;
    
    Int4 from = 0, to = query_length - lut_size;
    Int4 curr = from;
    Int4 scan_step = word_size - lut_size + 1;
    Uint4 index;
    BwtIntv ik, ok[4];
    Int4 ext_to = word_size - lut_size;
    Int4 ext_left, ext_right;
    Int4 i;
    Uint1 c;
    
    SearchInterval si;
    EBwtIntv ebi;
    
    std::vector<EBwtIntv> elist;
    
    Int8 total = 0;
    Uint8 scan_info[ext_to];
    
    for (curr = from; curr <= to; curr += scan_step)
    {
        index = GetSeedIndex(query + curr, lut_size);
        ik = ftable[index];
        if (ik.x[2] == 0) continue;
   
        i = curr - 1;
        ext_left = ext_right = 0;
        elist.clear();
        
        memset(scan_info, 0, sizeof(Uint8) * ext_to);
        
        while (ext_left < ext_to && i >= from)
        {
            c = query[i];
            Extend(&ik, ok, 1);
            if (ok[c].x[2] == 0) 
                break;
            if (ok[c].x[2] != ik.x[2])
            {
                ebi.intv = ik;
                ebi.ext_left = ext_left;
                elist.push_back(ebi);
            }

            ik = ok[c];
            ++ext_left;
            --i;
        }

		if (ik.x[2] != 0)
		{
			ebi.ext_left = ext_left;
			ebi.intv = ik;
			elist.push_back(ebi);
		}   
        
        std::vector<EBwtIntv>::reverse_iterator iter = elist.rbegin();
        ext_right = 0;
        i = curr + lut_size;
        ext_left = iter->ext_left + lut_size;
        ik = iter->intv;
        while (ext_left + ext_right < word_size && i < query_length)
        {
            c = 3 - query[i];
            Extend(&ik, ok, 0);
            if (ok[c].x[2] == 0) break;
            ik = ok[c];  
            scan_info[i - curr - lut_size] = ik.x[2];
            ++i;
            ++ext_right;
        }
        
        if (ext_left + ext_right == word_size)
        {
            si.query_id = context / 2;
            si.q_off = curr + lut_size - ext_left + q_start;
            si.k = ik.x[0];
            si.l = ik.x[2];   
            intvs.push_back(si);
        }
        
        ++iter;
        while (iter != elist.rend())
        {
            i = curr + lut_size;
            ext_left = iter->ext_left + lut_size;
            ik = iter->intv;
            ext_right = 0;
            while (ext_left + ext_right < word_size && i < query_length)
            {
                c = 3 - query[i];
                Extend(&ik, ok, 0);
                if (ok[c].x[2] == 0) break;
                ik = ok[c];
            
            if (ik.x[2] == scan_info[i - curr - lut_size]) break;
            scan_info[i - curr - lut_size] = ik.x[2];
            ++i;
            ++ext_right;
            }
            
            if (ext_left + ext_right == word_size)
            {
                si.query_id = context / 2;
                si.q_off = curr + lut_size - ext_left + q_start;
                si.k = ik.x[0];
                si.l = ik.x[2]; 
                intvs.push_back(si);  
            }
            ++iter;
        }
    }

    return 0;
}



















