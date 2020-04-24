#include "mem_finder.h"

#include "mem_finder.h"

#include "../corelib/ksort.h"
#include "hash_list_bucket_sort.h"

#define PV_BIT_SET(pv, pos) (((pv)[(pos)>>3]) |= (U8_ONE<<((pos)&7)))
#define PV_BIT_GET(pv, pos) (((pv)[(pos)>>3]) & (U8_ONE<<((pos)&7)))

#define kmif_hash_lt(a, b) ( \
    ((a).hash < (b).hash) \
    || \
    ((a).hash == (b).hash && (a).offset < (b).offset) \
)

KSORT_INIT(kmif_hash_lt, KmerInfo, kmif_hash_lt);

SmallLookupTable*
SmallLookupTableNew(int kmer_size)
{
    SmallLookupTable* slt = (SmallLookupTable*)calloc(1, sizeof(SmallLookupTable));
    const u64 kmer_cnt = 1<<(2*kmer_size);
    u64 pv_bytes = kmer_cnt / 8;
    slt->pv = (u8*)calloc(pv_bytes, 1);
    kv_init(slt->kmif_list);
    slt->kmer_stats = (int*)malloc(sizeof(int) * kmer_cnt);
    slt->kmer_size = kmer_size;
    slt->cv = (u8*)calloc(pv_bytes, 1);
    slt->seed_count = (int*)malloc(sizeof(int) * kmer_cnt);
    return slt;
}

SmallLookupTable*
SmallLookupTableFree(SmallLookupTable* slt)
{
    free(slt->pv);
    kv_destroy(slt->kmif_list);
    free(slt->kmer_stats);
    free(slt->cv);
    free(slt->seed_count);
    free(slt);
    return NULL;
}

MaximalExactMatchWorkData*
MaximalExactMatchWorkDataNew(int kmer_size, int window_size, int mem_size, int chain_score)
{
    MaximalExactMatchWorkData* data = 
        (MaximalExactMatchWorkData*)calloc(1, sizeof(MaximalExactMatchWorkData));
    kv_init(data->qry_kmif_list);
    kv_reserve(KmerInfo, data->qry_kmif_list, 100000);
    kv_init(data->ref_kmif_list);
    kv_reserve(KmerInfo, data->ref_kmif_list, 100000);
    kv_init(data->fwd_chain_seed_list);
    kv_reserve(ChainSeed, data->fwd_chain_seed_list, 100000);
    kv_init(data->rev_chain_seed_list);
    kv_reserve(ChainSeed, data->rev_chain_seed_list, 100000);
    kv_init(data->fwd_init_hit_list);
    kv_init(data->rev_init_hit_list);
    kv_init(data->init_hit_list);
    data->slt = SmallLookupTableNew(kmer_size);
    data->chain_data = ChainWorkDataNew(1, chain_score);
    data->kmer_size = kmer_size;
    data->window_size = window_size;
    data->mem_size = mem_size;
    data->fwd_ref = NULL;
    data->rev_ref = NULL;
    data->ref_size = 0;
    return data;
}

MaximalExactMatchWorkData*
MaximalExactMatchWorkDataFree(MaximalExactMatchWorkData* data)
{
    kv_destroy(data->qry_kmif_list);
    kv_destroy(data->ref_kmif_list);
    kv_destroy(data->fwd_chain_seed_list);
    kv_destroy(data->rev_chain_seed_list);
    kv_destroy(data->fwd_init_hit_list);
    kv_destroy(data->rev_init_hit_list);
    kv_destroy(data->init_hit_list);
    SmallLookupTableFree(data->slt);
    ChainWorkDataFree(data->chain_data);
    free(data);
    return NULL;
}

static int
build_subseq_kmif_list(const u8* sequence,
    const int seq_strand,
    const size_t seq_from,
    const size_t seq_to,
    const size_t seq_len, 
    const BOOL use_rc_seq,
    int kmer_size, 
    int window_size, 
    vec_kmif* kmif_list)
{
    KmerInfo kmif = { 0, 0, 0, seq_strand };
    if (seq_to - seq_from < kmer_size) return 0;
    const u8* s = sequence;
    const int intersect = kmer_size > window_size;
	u64 intersect_mask = 0;
	int stride = kmer_size - window_size;
	if (intersect) intersect_mask = (U64_ONE << (stride << 1)) - 1;

	if (!intersect) {
		for (size_t j = seq_from; j <= seq_to - kmer_size; j += window_size) {
			u64 hash = 0;
			for (int k = 0; k < kmer_size; ++k) {
                u8 c = use_rc_seq ? (3 - s[seq_len - 1 - j - k]) : s[j + k];
                hbn_assert(c >= 0 && c < 4, "c = %d,seq_from = %zu, seq_to = %zu, kmer_size = %d, window_size = %d, j = %zu, k = %d",
                    c, seq_from, seq_to, kmer_size, window_size, j, k);
				hash = (hash << 2) | c;
			}
            kmif.hash = hash;
            kmif.offset = j;
            kv_push(KmerInfo, *kmif_list, kmif);
		}
	} else {
		u64 hash = 0;
		for (int j = seq_from; j < seq_from + kmer_size; ++j) {
            u8 c = use_rc_seq ? (3 - s[seq_len - 1 - j]) : s[j];
            hbn_assert(c >= 0 && c < 4, "c = %d, seq_from = %zu, seq_to = %zu, kmer_size = %d, window_size = %d, j = %d",
                c, seq_from, seq_to, kmer_size, window_size, j);
			hash = (hash << 2) | c;
		}
		kmif.hash = hash;
        kmif.offset = seq_from;
        kv_push(KmerInfo, *kmif_list, kmif);
		for (u64 j = seq_from + window_size; j <= seq_to - kmer_size; j += window_size) {
			hash &= intersect_mask;
			for (int k = stride; k < kmer_size; ++k) {
                u8 c = use_rc_seq ? (3 - s[seq_len - 1 - j - k]) : s[j + k];
                hbn_assert(c >= 0 && c < 4, "c = %d,seq_from = %zu, seq_to = %zu, kmer_size = %d, window_size = %d, j = %zu, k = %d",
                    c, seq_from, seq_to, kmer_size, window_size, j, k);
				hash = (hash << 2) | c;
			}
			kmif.hash = hash;
            kmif.offset = j;
           // HBN_LOG("add %d", kmif.offset);
            kv_push(KmerInfo, *kmif_list, kmif);
		}
	}

    return kv_size(*kmif_list);
}


static int
build_kmif_list(const u8* sequence, 
    const int seq_strand,
    const size_t seq_len,
    const BOOL use_rc_seq,
    int kmer_size, 
    int window_size, 
    vec_kmif* kmif_list,
    const int SL,
    const int SR)
{
    int s = 0;
    const int n = seq_len;
    while (s < n) {
        int e = s + SL;
        e = hbn_min(e, n);
        int x = n - e;
        if (x <= SL + kmer_size * 2) e = n;
        build_subseq_kmif_list(sequence, seq_strand, s, e, seq_len, use_rc_seq, kmer_size, window_size, kmif_list);
        s = e + SR;
    }
    return kv_size(*kmif_list);
}

u64 kmer_info_list_hash_extractor(void* list, const u64 i)
{
    KmerInfo* kil = (KmerInfo*)list;
    return kil[i].hash;
}

u64 kmer_info_list_offset_extractor(void* list, const u64 i)
{
    KmerInfo* kil = (KmerInfo*)list;
    return kil[i].offset;    
}

void kmer_info_list_set_value(void* src, const u64 src_idx, void* dst, const u64 dst_idx)
{
    KmerInfo* srcl = (KmerInfo*)src;
    KmerInfo* dstl = (KmerInfo*)dst;
    dstl[dst_idx] = srcl[src_idx];
}

static void
sort_kmif_list(vec_kmif* kmif_list)
{
    size_t n = kv_size(*kmif_list);
    KmerInfo* a = kv_data(*kmif_list);
    if (n == 0) return;
    ks_introsort_kmif_hash_lt(n, a);

    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && a[j].hash == a[i].hash) ++j;
        a[i].occ = j - i;
        //HBN_LOG("offset = %d, occ = %d", a[i].offset, a[i].occ);
        i = j;
    }
}

void
validate_mem(HBN_LOG_PARAMS_GENERIC,
    const u8* read, 
    const u8* subject,
    const ChainSeed* cdpsa,
    const int cdpsc)
{
    return;
    //HBN_LOG("validating mems from [%s, %s, %d]", HBN_LOG_ARGS_GENERIC);
    for (int i = 0; i < cdpsc; ++i) {
        ChainSeed s = cdpsa[i];
        fprintf(stderr, "\tvalidating %d\t%d\t%d\t%d\t%d\n", i, s.qoff, s.soff, s.length, s.sdir);
        int qi = s.qoff;
        int si = s.soff;
        for (int k = 0; k < s.length; ++k, ++qi, ++si) {
            hbn_assert(read[qi] == subject[si], "[%s, %s, %d] at (%d, %d, %d): qi = %d, si = %d", 
            HBN_LOG_ARGS_GENERIC, s.qoff, s.soff, s.length, qi, si);
        }
    }
}

static int
left_extend(const u8* query,
    const int qsize,
    const u8* target,
    int qoff,
    int toff)
{
    int n = 0;
    while (qoff && toff) {
        --qoff;
        --toff;
        u8 qc = query[qoff];
        if (qc != target[toff]) break;
        ++n;
    }
    return n;
}

static int
right_extend(const u8* query,
    const u8* target,
    int qoff,
    int toff,
    int qsize,
    int tsize)
{
    int n = 0;
    while (qoff < qsize && toff < tsize) {
        u8 qc = query[qoff];
        if (qc != target[toff]) break;
        ++n;
        ++qoff;
        ++toff;
    }
    return n;
}

static void
extend_kmer_match(SmallLookupTable* slt,
    vec_chain_seed* mem_list,
    const u8* query,
    const int qsize,
    const u8* target,
    const int tsize,
    const int min_mem_size)
{
    size_t n = kv_size(*mem_list);
    ChainSeed* a = kv_data(*mem_list);
    if (n == 0) return;
    ks_introsort_chain_seed_soff_lt(n, a);
    const u8* q = query;
    const u8* t = target;
    int ql, qr, tl, tr;
    for (size_t i = 0; i < n; ++i) {
        if (a[i].length == 0) continue;
        ql = a[i].qoff;
        tl = a[i].soff;
        qr = ql + a[i].length;
        tr = tl + a[i].length;
        int e = left_extend(q, qsize, t, ql, tl);
        hbn_assert(ql >= e);
        ql -= e;
        hbn_assert(tl >= e);
        tl -= e;
        e = right_extend(q, t, qr, tr, qsize, tsize);
        hbn_assert(qr + e <= qsize);
        qr += e;
        hbn_assert(tr + e <= tsize);
        tr += e;
        a[i].qoff = ql;
        a[i].soff = tl;
        a[i].length = qr - ql;

        size_t j = i + 1;
        while (j < n && a[j].soff < tr) {
            if (a[j].qoff < qr && a[j].soff < tr
                &&
                qr - a[j].qoff == tr - a[j].soff) {
                a[j].length = 0;
            }
            ++j;
        }

        if (a[i].length < min_mem_size) a[i].length = 0;
    }

    size_t i = 0;
    for (size_t k = 0; k < n; ++k) {
        if (a[k].length) a[i++] = a[k];
    }
    kv_resize(ChainSeed, *mem_list, i);
}

void
SmallLookupTableBuild(SmallLookupTable* slt, 
    const u8* seq, 
    int seq_len)
{
    kv_clear(slt->kmif_list);
    build_kmif_list(seq, FWD, seq_len, FALSE, slt->kmer_size, 1, &slt->kmif_list, seq_len, 100);
    sort_kmif_list(&slt->kmif_list);
    int n = kv_size(slt->kmif_list);
    int i = 0;
    KmerInfo* kmif_array = kv_data(slt->kmif_list);
    const u64 max_kmer_hash = 1<<(2*slt->kmer_size);
    const u64 pv_bytes = max_kmer_hash / 8;
    memset(slt->pv, 0, pv_bytes);
    while (i < n) {
        if (kmif_array[i].occ<5) {
            u64 hash = kmif_array[i].hash;
		hbn_assert(hash < max_kmer_hash);
            PV_BIT_SET(slt->pv, hash);
            slt->kmer_stats[hash] = i;
        }
        i += kmif_array[i].occ;
    }
}

int
SmallLookupTableLookup(SmallLookupTable* slt, const u32 hash)
{
    if (!PV_BIT_GET(slt->pv, hash)) return -1;
    int p = slt->kmer_stats[hash];
    return p;
}

void 
MaximalExactMatchWorkData_Init(MaximalExactMatchWorkData* data, 
    const u8* fwd_ref,
    const u8* rev_ref, 
    const int ref_size)
{
    SmallLookupTableBuild(data->slt, fwd_ref, ref_size);
    kv_clear(data->fwd_chain_seed_list);
    kv_clear(data->rev_chain_seed_list);
    kv_clear(data->fwd_init_hit_list);
    kv_clear(data->rev_init_hit_list);
    kv_clear(data->init_hit_list);
    data->fwd_ref = fwd_ref;
    data->rev_ref = rev_ref;
    data->ref_size = ref_size;
}

void
collect_kmer_match(MaximalExactMatchWorkData* data, 
    ChainWorkData* chain_data,
    const u8* read, 
    const int read_size,
    const int read_strand)
{
    kv_clear(data->qry_kmif_list);
    if (read_strand == FWD || read_strand == F_R) {
        build_kmif_list(read, 
            FWD,
            read_size, 
            FALSE,
            data->kmer_size, 
            data->window_size, 
            &data->qry_kmif_list,
            read_size,
            0);
    }
    if (read_strand == REV || read_strand == F_R) {
        build_kmif_list(read, 
            REV,
            read_size, 
            TRUE,
            data->kmer_size, 
            data->window_size, 
            &data->qry_kmif_list,
            read_size,
            0);        
    }

    SmallLookupTable* slt = data->slt;
    kv_clear(chain_data->fwd_seeds);
    kv_clear(chain_data->rev_seeds);
    ChainSeed seed;
    seed.length = data->kmer_size;
    //HBN_LOG("number of query kmer: %zu", kv_size(data->qry_kmif_list));
	u64 max_kmer_hash = U64_ONE << (2*data->kmer_size);
    u64 cv_bytes = 1 << (2*data->kmer_size);
    cv_bytes /= 8;
    memset(slt->cv, 0, cv_bytes);
    for (size_t i = 0; i < kv_size(data->qry_kmif_list); ++i) {
        //seed.qoff = data->window_size * i;
        KmerInfo ki = kv_A(data->qry_kmif_list, i);
	hbn_assert(ki.hash < max_kmer_hash);
        int sdir;
        if (ki.strand == FWD) {
            sdir = FWD;
            seed.qoff = ki.offset;
        } else {
            sdir = REV;
            seed.qoff = read_size - data->kmer_size - ki.offset;
        }
        int k = SmallLookupTableLookup(slt, ki.hash);
        if (k == -1) continue;
        int cnt = kv_A(slt->kmif_list, k).occ;
        for (int j = 0; j < cnt; ++j) {
            KmerInfo ski = kv_A(slt->kmif_list, k+j);
		hbn_assert(ski.hash < max_kmer_hash);
            //HBN_LOG("ski offset = %d, occ = %d", ski.offset, ski.occ);
            if (ki.hash != ski.hash) {
                HBN_LOG("ref %s: i = %d, k = %d, j = %d, [%zu, %d, %d] v.s. [%zu, %d, %d], ref_size = %d, kmer_size = %d, window = %d, mem_size = %d",
                    data->ref_name, i, k, j, ki.hash, ki.occ, ki.occ, ski.hash, ski.occ, ski.offset, data->ref_size, data->kmer_size, data->window_size, data->mem_size);
            }
            hbn_assert(ki.hash == ski.hash);
            hbn_assert(ski.strand == FWD);
            if (sdir == FWD) {
                seed.soff = ski.offset;
                seed.hash = ki.hash;
                seed.sdir = FWD;
            } else {
                seed.soff = data->ref_size - data->kmer_size - ski.offset;
                seed.hash = ki.hash;
                seed.sdir = REV;
            }
            //HBN_LOG("find seed %d\t%d\t%d", seed.qoff, seed.soff, seed.length);
            if (sdir == FWD) {
                kv_push(ChainSeed, chain_data->fwd_seeds, seed);
            } else {
                kv_push(ChainSeed, chain_data->rev_seeds, seed);
            }
        }
    }
}

static void
remove_repetative_seeds(vec_chain_seed* seed_list)
{
    ChainSeed* seed_array = kv_data(*seed_list);
    int seed_count = kv_size(*seed_list);
    ks_introsort_chain_seed_soff_lt(seed_count, seed_array);
    int i = 0;
    while (i < seed_count) {
        int j = i + 1;
        while (j < seed_count && seed_array[i].soff == seed_array[j].soff) ++j;
        if (j - i > 5) {
            for (int k = i; k < j; ++k) seed_array[k].length = 0;
            //HBN_LOG("removing duplicated %d seeds", j - i);
        }
        i = j;
    }
    int k = 0;
    for (i = 0; i < seed_count; ++i) {
        if (seed_array[i].length > 0) seed_array[k++] = seed_array[i];
    }
    kv_resize(ChainSeed, *seed_list, k);
}

int
MaximalExactMatchWorkData_FindCandidates(
    MaximalExactMatchWorkData* data, 
    const u8* read, 
    const int read_size,
    const int read_strand)
{
    kv_clear(data->fwd_chain_seed_list);
    kv_clear(data->rev_chain_seed_list);
    kv_clear(data->fwd_init_hit_list);
    kv_clear(data->rev_init_hit_list);
    kv_clear(data->init_hit_list);
    
    const int verbose = 0;
    ChainWorkData* chain_data = data->chain_data;
    collect_kmer_match(data, chain_data, read, read_size, read_strand);

    extend_kmer_match(data->slt, 
        &chain_data->fwd_seeds, 
        read, 
        read_size, 
        data->fwd_ref, 
        data->ref_size, 
        data->mem_size);

    extend_kmer_match(data->slt,
        &chain_data->rev_seeds,
        read,
        read_size,
        data->rev_ref,
        data->ref_size,
        data->mem_size);

    remove_repetative_seeds(&chain_data->fwd_seeds);
    remove_repetative_seeds(&chain_data->rev_seeds);

    chaining_find_candidates(chain_data, 
        kv_data(chain_data->fwd_seeds),
        kv_size(chain_data->fwd_seeds),
        FWD,
        &data->fwd_init_hit_list, 
        &data->fwd_chain_seed_list);

#if 0
HBN_LOG("fwd mems:");
for (size_t i = 0; i < kv_size(chain_data->fwd_seeds); ++i) {
    ChainSeed s = kv_A(chain_data->fwd_seeds, i);
    fprintf(stderr, "\t%d\t%d\t%d\t%d\n", i, s.qoff, s.soff, s.length);
}
#endif

    chaining_find_candidates(chain_data, 
        kv_data(chain_data->rev_seeds),
        kv_size(chain_data->rev_seeds),
        REV,
        &data->rev_init_hit_list, 
        &data->rev_chain_seed_list);

#if 0
HBN_LOG("rev mems:");
for (size_t i = 0; i < kv_size(chain_data->rev_seeds); ++i) {
    ChainSeed s = kv_A(chain_data->rev_seeds, i);
    fprintf(stderr, "\t%d\t%d\t%d\t%d\n", i, s.qoff, s.soff, s.length);
}
#endif

    kv_push_v(HbnInitHit, 
        data->init_hit_list, 
        kv_data(data->fwd_init_hit_list), 
        kv_size(data->fwd_init_hit_list));
    kv_push_v(HbnInitHit,
        data->init_hit_list,
        kv_data(data->rev_init_hit_list),
        kv_size(data->rev_init_hit_list));

    ks_introsort_init_hit_score_gt(kv_size(data->init_hit_list), kv_data(data->init_hit_list));
    return !kv_empty(data->init_hit_list);
}
