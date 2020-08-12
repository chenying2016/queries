#include "primer_map_hit_finder.h"

#include "../../corelib/cstr_util.h"
#include "../../corelib/ksort.h"
#include "../../ncbi_blast/c_ncbi_blast_aux.h"

#include "sort_primer_map_seeds.h"

#define pm_word_hash_lt(a, b) ((a).hash < (b).hash)
KSORT_INIT(pm_word_hash_lt, PrimerMapWord, pm_word_hash_lt);

#define pm_word_key_lt(a, b) ((a).key < (b).key)
KSORT_INIT(pm_word_key_lt, PrimerMapWord, pm_word_key_lt);

#define pm_seed_key_lt(a, b) ((a).key < (b).key)
KSORT_INIT(pm_seed_key_lt, PrimerMapSeed, pm_seed_key_lt);

static const int kMaxKmerOcc = 200;

PrimerMapHitFindData*
PrimerMapHitFindDataNew(int word_size, int word_stride, int chain_score)
{
    PrimerMapHitFindData* data = (PrimerMapHitFindData*)calloc(1, sizeof(PrimerMapHitFindData));
    data->word_size = word_size;
    data->word_stride = word_stride;
    data->chain = ChainWorkDataNew_PrimerMap(1, chain_score);
    data->num_primer_contexts = 0;
    kv_init(data->primer_context_list);
    kv_init(data->primer_word_list);
    data->num_query_contexts = 0;
    kv_init(data->query_context_list);
    kv_init(data->query_word_list);
    kv_init(data->seed_list);
    kv_init(data->hit_list);
    kv_init(data->hit_seed_list);
    return data;
}

PrimerMapHitFindData*
PrimerMapHitFindDataFree(PrimerMapHitFindData* data)
{
    if (!data) return NULL;
    data->chain = ChainWorkDataFree(data->chain);
    kv_destroy(data->primer_context_list);
    kv_destroy(data->primer_word_list);
    kv_destroy(data->query_context_list);
    kv_destroy(data->query_word_list);
    kv_destroy(data->seed_list);
    kv_destroy(data->hit_list);
    kv_destroy(data->hit_seed_list);
    sfree(data);
    return NULL;
}

//////////////////////

const u8*
PrimerMapHitFindData_ExtractPrimer(PrimerMapHitFindData* data, int context)
{
    hbn_assert(context < data->num_primer_contexts);
    return kv_A(data->primer_context_list, context).sequence;
}

int
PrimerMapHitFindData_PrimerSize(PrimerMapHitFindData* data, int context)
{
    hbn_assert(context < data->num_primer_contexts);
    return kv_A(data->primer_context_list, context).size;
}

const char*
PrimerMapHitFindData_PrimerName(PrimerMapHitFindData* data, int context)
{
    hbn_assert(context < data->num_primer_contexts);
    return kv_A(data->primer_context_list, context).name;    
}

void
PrimerMapHitFindData_CleanPrimerData(PrimerMapHitFindData* data)
{
    data->num_primer_contexts = 0;
    kv_clear(data->primer_context_list);
    kv_clear(data->primer_word_list);
}

void
PrimerMapHitFindData_AddOnePrimer(PrimerMapHitFindData* data,
    const int primer_index,
    const u8* primer,
    const int primer_size,
    const char* primer_name)
{
    PrimerMapContextInfo ctx;
    ctx.context = data->num_primer_contexts;
    ctx.seq_index = primer_index;
    ctx.name = primer_name;
    ctx.size = primer_size;
    ctx.sequence = primer;
    kv_push(PrimerMapContextInfo, data->primer_context_list, ctx);
    data->num_primer_contexts++;
}

////////////////////////////

const u8*
PrimerMapHitFindData_ExtractQuery(PrimerMapHitFindData* data, int context)
{
    hbn_assert(context < data->num_query_contexts);
    return kv_A(data->query_context_list, context).sequence;
}

int
PrimerMapHitFindData_QuerySize(PrimerMapHitFindData* data, int context)
{
    hbn_assert(context < data->num_query_contexts);
    return kv_A(data->query_context_list, context).size;    
}

const char*
PrimerMapHitFindData_QueryName(PrimerMapHitFindData* data, int context)
{
    hbn_assert(context < data->num_query_contexts);
    return kv_A(data->query_context_list, context).name;    
}

void
PrimerMapHitFindData_CleanQueryData(PrimerMapHitFindData* data)
{
    data->num_query_contexts = 0;
    kv_clear(data->query_context_list);
    kv_clear(data->query_word_list);
}

void
PrimerMapHitFindData_AddOneQuery(PrimerMapHitFindData* data,
    const int query_index,
    const u8* query,
    const int query_size,
    const char* query_name)
{
    PrimerMapContextInfo ctx;
    ctx.seq_index = query_index;
    ctx.context = data->num_query_contexts;
    ctx.name = query_name;
    ctx.size = query_size;
    ctx.sequence = query;
    kv_push(PrimerMapContextInfo, data->query_context_list, ctx);
    data->num_query_contexts++;
}

//////////////////////

static int
build_subseq_word_list(const u8* sequence,
    const int context,
    const size_t seq_from,
    const size_t seq_to,
    const size_t seq_len,
    int kmer_size, 
    int window_size, 
    vec_pm_word* word_list)
{
    PrimerMapWord word;
    word.context = context;
    word.hash = 0;
    word.key = 0;
    word.offset = 0;
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
                u8 c = s[j + k];
                hbn_assert(c >= 0 && c < 4, "c = %d,seq_from = %zu, seq_to = %zu, kmer_size = %d, window_size = %d, j = %zu, k = %d",
                    c, seq_from, seq_to, kmer_size, window_size, j, k);
				hash = (hash << 2) | c;
			}
            word.hash = hash;
            word.offset = j;
            kv_push(PrimerMapWord, *word_list, word);
		}
	} else {
		u64 hash = 0;
		for (int j = seq_from; j < seq_from + kmer_size; ++j) {
            u8 c = s[j];
            hbn_assert(c >= 0 && c < 4, "c = %d, seq_from = %zu, seq_to = %zu, kmer_size = %d, window_size = %d, j = %d",
                c, seq_from, seq_to, kmer_size, window_size, j);
			hash = (hash << 2) | c;
		}
		word.hash = hash;
        word.offset = seq_from;
        kv_push(PrimerMapWord, *word_list, word);
		for (u64 j = seq_from + window_size; j <= seq_to - kmer_size; j += window_size) {
			hash &= intersect_mask;
			for (int k = stride; k < kmer_size; ++k) {
                u8 c = s[j + k];
                hbn_assert(c >= 0 && c < 4, "c = %d,seq_from = %zu, seq_to = %zu, kmer_size = %d, window_size = %d, j = %zu, k = %d",
                    c, seq_from, seq_to, kmer_size, window_size, j, k);
				hash = (hash << 2) | c;
			}
			word.hash = hash;
            word.offset = j;
            kv_push(PrimerMapWord, *word_list, word);
		}
	}

    return kv_size(*word_list);
}

static int
build_word_list(const u8* sequence, 
    const int context,
    const size_t seq_len,
    int kmer_size, 
    int window_size, 
    vec_pm_word* word_list,
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
        build_subseq_word_list(sequence, context, s, e, seq_len, kmer_size, window_size, word_list);
        s = e + SR;
    }
    return kv_size(*word_list);
}

static int
s_proceed_to_next_word_hash_idx(PrimerMapWord* a, int n, int i)
{
    if (i >= n) return i;
    int hash = a[i].hash;
    ++i;
    while (i < n) {
        if (a[i].hash != hash) break;
        ++i;
    }
    return i;
}

static int
s_proceed_to_next_word_key_idx(PrimerMapWord* a, int n, int i)
{
    if (i >= n) return i;
    u64 key = a[i].key;
    ++i;
    while (i < n) {
        if (a[i].key != key) break;
        ++i;
    }
    return i;
}

void
PrimerMapHitFindData_BuildPrimerWordList(PrimerMapHitFindData* data)
{
    for (int i = 0; i < data->num_primer_contexts; ++i) {
        const u8* primer = PrimerMapHitFindData_ExtractPrimer(data, i);
        const int primer_size = PrimerMapHitFindData_PrimerSize(data, i);
        //HBN_LOG("primer %d size: %d", i, primer_size);
        build_word_list(primer, i, primer_size,
            data->word_size, 1, &data->primer_word_list, primer_size, 0);
    }
    PrimerMapWord* wa = kv_data(data->primer_word_list);
    int wc = kv_size(data->primer_word_list);
    ks_introsort_pm_word_hash_lt(wc, wa);
    for (int i = 0; i < wc - 1; ++i) {
        PrimerMapWord w1 = wa[i];
        PrimerMapWord w2 = wa[i+1];
        hbn_assert(w1.hash <= w2.hash);
    }
}

void
PrimerMapHitFindData_BuildQueryWordList(PrimerMapHitFindData* data)
{
    for (int i = 0; i < data->num_query_contexts; ++i) {
        const u8* query = PrimerMapHitFindData_ExtractQuery(data, i);
        const int query_size = PrimerMapHitFindData_QuerySize(data, i);
        build_word_list(query, i, query_size,
            data->word_size, data->word_stride, &data->query_word_list, query_size, 0);
    }    
    PrimerMapWord* wa = kv_data(data->query_word_list);
    int wc = kv_size(data->query_word_list);
    for (int i = 0; i < wc; ++i) {
        u64 hash = wa[i].hash;
        u64 context = wa[i].context;
        u64 key = (hash << 32) | context;
        wa[i].key = key;
    }
    ks_introsort_pm_word_key_lt(wc, wa);
    for (int i = 0; i < wc - 1; ++i) {
        PrimerMapWord w1 = wa[i];
        PrimerMapWord w2 = wa[i+1];
        hbn_assert((w1.hash < w2.hash) || (w1.hash == w2.hash && w1.context <= w2.context));
    }
}

static void
s_collect_seeds(PrimerMapHitFindData* data)
{
    PrimerMapWord* p_wa = kv_data(data->primer_word_list);
    const int p_wc = kv_size(data->primer_word_list);
    PrimerMapWord* q_wa = kv_data(data->query_word_list);
    const int q_wc = kv_size(data->query_word_list);
    int p_wi = 0;
    int q_wi = 0;
    PrimerMapSeed seed;
    seed.key = 0;
    seed.length = data->word_size;
    //HBN_LOG("primer word: %d, query word: %d", pwc, swc);

    while (q_wi < q_wc && p_wi < p_wc) {
        while (q_wi < q_wc && q_wa[q_wi].hash < p_wa[p_wi].hash) {
            q_wi = s_proceed_to_next_word_hash_idx(q_wa, q_wc, q_wi);
        }
        if (q_wi >= q_wc) break;
        while (p_wi < p_wc && p_wa[p_wi].hash < q_wa[q_wi].hash) {
            p_wi = s_proceed_to_next_word_hash_idx(p_wa, p_wc, p_wi);
        }
        if (p_wi >= p_wc) break;
        if (q_wa[q_wi].hash != p_wa[p_wi].hash) continue;
        int next_p_wi = s_proceed_to_next_word_hash_idx(p_wa, p_wc, p_wi);
        int next_q_wi = s_proceed_to_next_word_hash_idx(q_wa, q_wc, q_wi);

        int q_i = q_wi;
        while (q_i < next_q_wi) {
            int q_i_to = s_proceed_to_next_word_key_idx(q_wa, next_q_wi, q_i);
            if (q_i_to - q_i > kMaxKmerOcc) {
                q_i = q_i_to;
                continue;
            }
            for (int p_idx = p_wi; p_idx < next_p_wi; ++p_idx) {
                seed.primer_context = p_wa[p_idx].context;
                seed.primer_offset = p_wa[p_idx].offset;
                for (int q_idx = q_i; q_idx < q_i_to; ++q_idx) {
                    seed.query_context = q_wa[q_idx].context;
                    seed.query_offset = q_wa[q_idx].offset;
                    if (0 && seed.primer_offset == 0) {
                        HBN_LOG("**** qid = %d, qoff = %d, sid = %d, soff = %d", 
                            seed.primer_context, seed.primer_offset, 
                            seed.query_context, seed.query_offset);
                    }
                    kv_push(PrimerMapSeed, data->seed_list, seed);
                }
            }
            q_i = q_i_to;
        }
        q_wi = next_q_wi;
        p_wi = next_p_wi;
    }
    //HBN_LOG("number of kmer match: %zu", kv_size(data->seed_list));
}

static void
s_fill_init_hit_list_info(
    int primer_context,
    int primer_size,
    int query_context,
    int query_size,
    HbnInitHit* hit_array,
    int hit_count)
{
    for (int i = 0; i < hit_count; ++i) {
        HbnInitHit* hit = hit_array + i;
        hit->qid = primer_context;
        hit->qdir = FWD;
        hit->qsize = primer_size;
        hit->sid = query_context;
        hit->sdir = FWD;
        hit->ssize = query_size;
        ChainSeed* fcs = hit->chain_seed_array;
        ChainSeed* lcs = fcs + hit->chain_seed_count - 1;
        hit->qbeg = fcs->qoff;
        hit->sbeg = fcs->soff;
        hit->qend = lcs->qoff + lcs->length;
        hit->send = lcs->soff + lcs->length;
        hit->qoff = hit->qend;
        hit->soff = hit->send;
    }    
}

static void
s_reverse_init_hit_info(const HbnInitHit* src, HbnInitHit* dst)
{
    dst->qid = src->sid;
    dst->qdir = src->sdir;
    dst->qsize = src->ssize;
    dst->qbeg = src->sbeg;
    dst->qend = src->send;
    dst->qoff = src->soff;

    dst->sid = src->qid;
    dst->sdir = src->qdir;
    dst->ssize = src->qsize;
    dst->sbeg = src->qbeg;
    dst->send = src->qend;
    dst->soff = src->qoff;

    dst->score = src->score;
    dst->chain_seed_array = NULL;
    dst->chain_seed_count = src->chain_seed_count;
    dst->chain_seed_offset = src->chain_seed_offset;
}

static void
s_find_candidates_for_one_context(PrimerMapHitFindData* data, PrimerMapSeed* sa, int sc)
{
    const int primer_context = sa[0].primer_context;
    const int primer_size = PrimerMapHitFindData_PrimerSize(data, primer_context);
    const int query_context = sa[0].query_context;
    const int query_size = PrimerMapHitFindData_QuerySize(data, query_context);
    
    kv_clear(data->chain->fwd_seeds);
    ChainSeed cs;
    for (int i = 0; i < sc; ++i) {
        cs.hash = 0;
        cs.qoff = sa[i].primer_offset;
        cs.sdir = FWD;
        cs.soff = sa[i].query_offset;
        cs.length = sa[i].length;
        //HBN_LOG("%d\tqoff = %d, soff = %d", i, cs.qoff, cs.soff);
        kv_push(ChainSeed, data->chain->fwd_seeds, cs);
    }
    sort_chain_seed_soff_lt(kv_size(data->chain->fwd_seeds), kv_data(data->chain->fwd_seeds));

    kv_dinit(vec_init_hit, hit_list);
    kv_dinit(vec_chain_seed, hit_seed_list);
    chaining_find_candidates_primer_map(data->chain,
        kv_data(data->chain->fwd_seeds),
        kv_size(data->chain->fwd_seeds),
        FALSE,
        0,
        &hit_list,
        &hit_seed_list);  

    //HBN_LOG("number of hits: %zu", kv_size(hit_list));

    s_fill_init_hit_list_info(
        primer_context,
        primer_size,
        query_context,
        query_size,
        kv_data(hit_list),
        kv_size(hit_list));

    HbnInitHit new_hit;
    for (int i = 0; i < kv_size(hit_list); ++i) {
        HbnInitHit* hit = &kv_A(hit_list, i);
        s_reverse_init_hit_info(hit, &new_hit);
        new_hit.chain_seed_offset = kv_size(data->hit_seed_list);
        kv_push(HbnInitHit, data->hit_list, new_hit);

        ChainSeed* csa = (ChainSeed*)hit->chain_seed_array;
        int csc = hit->chain_seed_count;
        for (int j = 0; j < csc; ++j) {
            cs.hash = csa[j].hash;
            cs.length = csa[j].length;
            cs.qoff = csa[j].soff;
            cs.soff = csa[j].qoff;
            cs.sdir = FWD;
            kv_push(ChainSeed, data->hit_seed_list, cs);
        }
    }
    kv_destroy(hit_list);
    kv_destroy(hit_seed_list);
}

static void
s_find_candidates(PrimerMapHitFindData* data)
{
    PrimerMapSeed* seed_array = kv_data(data->seed_list);
    int seed_count = kv_size(data->seed_list);
    for (int i = 0; i < seed_count; ++i) {
        hbn_assert(seed_array[i].length > 0);
        u64 pid = seed_array[i].primer_context;
        u64 qid = seed_array[i].query_context;
        seed_array[i].key = (qid<<32) | pid;
    }
    ks_introsort_pm_seed_key_lt(seed_count, seed_array);
    for (int i = 0; i < seed_count - 1; ++i) {
        PrimerMapSeed c1 = seed_array[i];
        PrimerMapSeed c2 = seed_array[i+1];
        hbn_assert((c1.query_context < c2.query_context)
                    ||
                    (c1.query_context == c2.query_context && c1.primer_context <= c2.primer_context));
    }

    int i = 0;
    while (i < seed_count) {
        int j = i + 1;
        while (j < seed_count 
               && 
               seed_array[i].query_context == seed_array[j].query_context 
               && 
               seed_array[i].primer_context == seed_array[j].primer_context) {
            ++j;
        }
        s_find_candidates_for_one_context(data, seed_array + i, j - i);
        i = j;
    }
}

void
PrimerMapHitFindData_FindHits(PrimerMapHitFindData* data)
{
    kv_clear(data->seed_list);
    kv_clear(data->hit_list);
    kv_clear(data->hit_seed_list);

    s_collect_seeds(data);
    s_find_candidates(data);

    for (size_t i = 0; i < kv_size(data->hit_list); ++i) {
        kv_A(data->hit_list, i).chain_seed_array = kv_data(data->hit_seed_list) 
                                                   + 
                                                   kv_A(data->hit_list, i).chain_seed_offset;
    }
}