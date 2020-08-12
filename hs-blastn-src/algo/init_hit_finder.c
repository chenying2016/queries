#include "init_hit_finder.h"

#include "../corelib/ksort.h"

#define ihf_word_hash_strand_lt(a, b) (\
    ((a).hash < (b).hash) \
    || \
    ((a).hash == (b).hash && (a).strand < (b).strand) \
)
KSORT_INIT(ihf_word_hash_strand_lt, InitHitFindWord, ihf_word_hash_strand_lt);

#define ihf_seed_off_lt(a, b) ( \
    ((a).subject_offset < (b).subject_offset) \
    || \
    ((a).subject_offset == (b).subject_offset && (a).query_offset < (b).query_offset) \
)
KSORT_INIT(ihf_seed_off_lt, InitHitFindSeed, ihf_seed_off_lt);

#define chain_seed_off_lt(a, b) ( \
    ((a).soff < (b).soff) \
    || \
    ((a).soff == (b).soff && (a).qoff < (b).qoff) \
)
KSORT_INIT(chain_seed_off_lt, ChainSeed, chain_seed_off_lt);

static const int kMaxWordOcc = 8;
static const int kMaxSeedOcc = 8;

/////////////////

static int
build_subseq_word_list(const u8* sequence,
    const int seq_strand,
    const size_t seq_from,
    const size_t seq_to,
    const size_t seq_len, 
    const BOOL use_rc_seq,
    int kmer_size, 
    int window_size, 
    vec_ihf_word* word_list)
{
    InitHitFindWord word;
    word.hash = 0;
    word.offset = 0;
    word.strand = seq_strand;
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
            word.hash = hash;
            word.offset = j;
            kv_push(InitHitFindWord, *word_list, word);
		}
	} else {
		u64 hash = 0;
		for (int j = seq_from; j < seq_from + kmer_size; ++j) {
            u8 c = use_rc_seq ? (3 - s[seq_len - 1 - j]) : s[j];
            hbn_assert(c >= 0 && c < 4, "c = %d, seq_from = %zu, seq_to = %zu, kmer_size = %d, window_size = %d, j = %d",
                c, seq_from, seq_to, kmer_size, window_size, j);
			hash = (hash << 2) | c;
		}
		word.hash = hash;
        word.offset = seq_from;
        kv_push(InitHitFindWord, *word_list, word);
		for (u64 j = seq_from + window_size; j <= seq_to - kmer_size; j += window_size) {
			hash &= intersect_mask;
			for (int k = stride; k < kmer_size; ++k) {
                u8 c = use_rc_seq ? (3 - s[seq_len - 1 - j - k]) : s[j + k];
                hbn_assert(c >= 0 && c < 4, "c = %d,seq_from = %zu, seq_to = %zu, kmer_size = %d, window_size = %d, j = %zu, k = %d",
                    c, seq_from, seq_to, kmer_size, window_size, j, k);
				hash = (hash << 2) | c;
			}
			word.hash = hash;
            word.offset = j;
           // HBN_LOG("add %d", kmif.offset);
            kv_push(InitHitFindWord, *word_list, word);
		}
	}

    return kv_size(*word_list);
}

static int
build_word_list(const u8* sequence, 
    const int seq_strand,
    const size_t seq_len,
    const BOOL use_rc_seq,
    int kmer_size, 
    int window_size, 
    vec_ihf_word* word_list,
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
        build_subseq_word_list(sequence, seq_strand, s, e, seq_len, use_rc_seq, 
            kmer_size, window_size, word_list);
        s = e + SR;
    }
    return kv_size(*word_list);
}

static int
s_proceed_to_next_word_hash_idx(InitHitFindWord* a, int n, int i)
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

static void
s_collect_seeds(InitHitFindData* data, int subject_dir)
{
    InitHitFindWord* q_wa = kv_data(data->query_word_list);
    const int q_wc = kv_size(data->query_word_list);
    int q_wi = 0;
    InitHitFindWord* s_wa = kv_data(data->subject_word_list);
    const int s_wc = kv_size(data->subject_word_list);
    int s_wi = 0;
    kv_clear(data->fwd_subject_seed_list);
    kv_clear(data->rev_subject_seed_list);
    InitHitFindSeed seed;

    while (s_wi < s_wc && q_wi < q_wc) {
        while (s_wi < s_wc && s_wa[s_wi].hash < q_wa[q_wi].hash) {
            s_wi = s_proceed_to_next_word_hash_idx(s_wa, s_wc, s_wi);
        }
        if (s_wi >= s_wc) break;

        while (q_wi < q_wc && q_wa[q_wi].hash < s_wa[s_wi].hash) {
            q_wi = s_proceed_to_next_word_hash_idx(q_wa, q_wc, q_wi);
        }
        if (q_wi >= q_wc) break;

        if (q_wa[q_wi].hash != s_wa[s_wi].hash) continue;

        int next_q_wi = s_proceed_to_next_word_hash_idx(q_wa, q_wc, q_wi);
        int next_s_wi = s_proceed_to_next_word_hash_idx(s_wa, s_wc, s_wi);
        if (next_q_wi - q_wi > kMaxWordOcc) {
            q_wi = next_q_wi;
            s_wi = next_s_wi;
            continue;
        }

        int s_wi_d = s_wi;
        while (s_wi_d < next_s_wi && s_wa[s_wi_d].strand == FWD) ++s_wi_d;
        int r;
        r = (subject_dir == FWD || subject_dir == F_R)
            &&
            (s_wi_d - s_wi <= kMaxWordOcc)
            &&
            ((s_wi_d - s_wi) * (next_q_wi - q_wi) <= kMaxSeedOcc);
        if (r) {
            for (int qi = q_wi; qi < next_q_wi; ++qi) {
                seed.query_offset = q_wa[qi].offset;
                for (int si = s_wi; si < s_wi_d; ++si) {
                    hbn_assert(q_wa[qi].hash == s_wa[si].hash);
                    seed.subject_offset = s_wa[si].offset;
                    kv_push(InitHitFindSeed, data->fwd_subject_seed_list, seed);
                }
            }
        }
        r = (subject_dir == REV || subject_dir == F_R)
            &&
            (next_s_wi - s_wi_d <= kMaxWordOcc)
            &&
            ((next_s_wi - s_wi_d) * (next_q_wi - q_wi) <= kMaxSeedOcc);
        if (r) {
            for (int qi = q_wi; qi < next_q_wi; ++qi) {
                seed.query_offset = q_wa[qi].offset;
                for (int si = s_wi_d; si < next_s_wi; ++si) {
                    hbn_assert(q_wa[qi].hash == s_wa[si].hash);
                    seed.subject_offset = s_wa[si].offset;
                    kv_push(InitHitFindSeed, data->rev_subject_seed_list, seed);
                }
            }
        }

        q_wi = next_q_wi;
        s_wi = next_s_wi;
    }
}

static void
s_fill_init_hit_info(int qid, int qsize, int sid, int ssize, HbnInitHit* hit_array, int hit_count)
{
    for (int i = 0; i < hit_count; ++i) {
        HbnInitHit* hit = hit_array + i;
        hit->qid = qid;
        hit->qdir = FWD;
        hit->qsize = qsize;
        hit->sid = sid;
        hit->ssize = ssize;
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

void
print_init_hit_range(HbnInitHit* hit)
{
    ChainSeed* csa = (ChainSeed*)(hit->chain_seed_array);
    int csc = hit->chain_seed_count;
    int qoff = csa[0].qoff;
    int soff = csa[0].soff;
    int qend = csa[csc-1].qoff + csa[csc-1].length;
    int send = csa[csc-1].soff + csa[csc-1].length;
    fprintf(stderr, "[%d, %d, %d, %d, %d] x [%d, %d, %d, %d, %zu], score = %d\n",
        hit->qid, hit->qdir, qoff, qend, hit->qsize, 
        hit->sid, hit->sdir, soff, send, hit->ssize,
        hit->score);
}

static void
s_extract_mem(const u8* query,
    const u8* subject,
    const int qoff,
    const int soff,
    const int length,
    const int mem_size,
    vec_chain_seed* seed_list)
{
    //HBN_LOG("**** qoff = %d, soff = %d, length = %d", qoff, soff, length);
    int i = 0;
    int qi = qoff, si = soff;
    while (i < length) {
        while (i < length) {
            if (query[qi] == subject[si]) break;
            ++qi;
            ++si;
            ++i;
        }
        if (length - i < mem_size) break;

        hbn_assert(query[qi] == subject[si]);
        int j = i;
        while (j < length) {
            if (query[qi] != subject[si]) break;
            ++qi;
            ++si;
            ++j;
        }

        if (j - i >= mem_size) {
            ChainSeed seed;
            seed.hash = 0;
            seed.length = j - i;
            seed.qoff = qi - seed.length;
            seed.soff = si - seed.length;
            seed.sdir = FWD;
            //HBN_LOG("qoff = %d, soff = %d, length = %d", seed.qoff, seed.soff, seed.length);
            kv_push(ChainSeed, *seed_list, seed);
        }
        i = j;
    }
}

////////////////////////////////////////

InitHitFindData*
InitHitFindDataNew(int word_size, int word_stride, int chain_score)
{
    InitHitFindData* data = (InitHitFindData*)calloc(1, sizeof(InitHitFindData));
    data->word_size = word_size;
    data->word_stride = word_stride;

    kv_init(data->query_word_list);
    kv_init(data->subject_word_list);
    kv_init(data->fwd_subject_seed_list);
    kv_init(data->rev_subject_seed_list);

    data->chain = ChainWorkDataNew(1, chain_score);
    kv_init(data->hit_list);
    kv_init(data->hit_seed_list);

    return data;
}

InitHitFindData*
InitHitFindDataFree(InitHitFindData* data)
{
    kv_destroy(data->query_word_list);
    kv_destroy(data->subject_word_list);
    kv_destroy(data->fwd_subject_seed_list);
    kv_destroy(data->rev_subject_seed_list);
    data->chain = ChainWorkDataFree(data->chain);
    kv_destroy(data->hit_list);
    kv_destroy(data->hit_seed_list);
    free(data);
    return NULL;
}

void
InitHitFindData_AddQuery(InitHitFindData* data,
    const int query_oid,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size)
{
    hbn_assert(fwd_query != NULL);

    data->query_oid = query_oid;
    data->query_name = query_name;
    data->fwd_query = fwd_query;
    data->rev_query = rev_query;
    data->query_size = query_size;

    kv_clear(data->query_word_list);
    if (fwd_query) build_word_list(fwd_query, FWD, query_size, FALSE, 
        data->word_size, data->word_stride,
        &data->query_word_list, query_size, 0);

    InitHitFindWord* wa = kv_data(data->query_word_list);
    size_t wc = kv_size(data->query_word_list);
    ks_introsort_ihf_word_hash_strand_lt(wc, wa);
}

void
InitHitFindData_Init(InitHitFindData* data,
    const int subject_oid,
    const char* subject_name,
    const u8* fwd_subject,
    const u8* rev_subject,
    const int subject_size)
{
    data->subject_oid = subject_oid;
    data->subject_name = subject_name;
    data->fwd_subject = fwd_subject;
    data->rev_subject = rev_subject;
    data->subject_size = subject_size;

    kv_clear(data->subject_word_list);
    if (fwd_subject) build_word_list(fwd_subject, FWD, subject_size, FALSE,
        data->word_size, 1,
        &data->subject_word_list, subject_size, 0);
    if (rev_subject) build_word_list(rev_subject, REV, subject_size, FALSE,
        data->word_size, 1,
        &data->subject_word_list, subject_size, 0);

    InitHitFindWord* wa = kv_data(data->subject_word_list);
    size_t wc = kv_size(data->subject_word_list);
    ks_introsort_ihf_word_hash_strand_lt(wc, wa);
}

void
InitHitFindData_FindHits(InitHitFindData* data, int subject_dir)
{
    s_collect_seeds(data, subject_dir);
    kv_clear(data->hit_list);
    kv_clear(data->hit_seed_list);

    InitHitFindSeed* fsa = kv_data(data->fwd_subject_seed_list);
    int fsc = kv_size(data->fwd_subject_seed_list);
    kv_clear(data->chain->fwd_seeds);
    for (int i = 0; i < fsc; ++i) {
        ChainSeed cs;
        cs.hash = 0;
        cs.length = data->word_size;
        cs.qoff = fsa[i].query_offset;
        cs.soff = fsa[i].subject_offset;
        cs.sdir = FWD;
        kv_push(ChainSeed, data->chain->fwd_seeds, cs);
    }
    //HBN_LOG("fwd seeds: %zu", kv_size(data->chain->fwd_seeds));
    ks_introsort_chain_seed_off_lt(kv_size(data->chain->fwd_seeds), kv_data(data->chain->fwd_seeds));
    kv_dinit(vec_init_hit, hit_list);
    kv_dinit(vec_chain_seed, hit_seed_list);
    chaining_find_candidates(data->chain,
        kv_data(data->chain->fwd_seeds),
        kv_size(data->chain->fwd_seeds),
        FALSE,
        FWD,
        &hit_list,
        &hit_seed_list); 
    s_fill_init_hit_info(data->query_oid, data->query_size, 
        data->subject_oid, data->subject_size, 
        kv_data(hit_list), kv_size(hit_list));
    for (size_t i = 0; i < kv_size(hit_list); ++i) {
        HbnInitHit* hit = &kv_A(hit_list, i);
        hit->chain_seed_offset = kv_size(data->hit_seed_list);
        kv_push(HbnInitHit, data->hit_list, *hit);
        kv_push_v(ChainSeed, data->hit_seed_list, hit->chain_seed_array, hit->chain_seed_count);
    }

    InitHitFindSeed* rsa = kv_data(data->rev_subject_seed_list);
    int rsc = kv_size(data->rev_subject_seed_list);
    kv_clear(data->chain->rev_seeds);
    for (int i = 0; i < rsc; ++i) {
        ChainSeed cs;
        cs.hash = 0;
        cs.length = data->word_size;
        cs.qoff = rsa[i].query_offset;
        cs.soff = rsa[i].subject_offset;
        cs.sdir = REV;
        kv_push(ChainSeed, data->chain->rev_seeds, cs);
    }
    //HBN_LOG("rev seeds: %zu", kv_size(data->chain->rev_seeds));
    ks_introsort_chain_seed_off_lt(kv_size(data->chain->rev_seeds), kv_data(data->chain->rev_seeds));
    kv_clear(hit_list);
    kv_clear(hit_seed_list);
    chaining_find_candidates(data->chain,
        kv_data(data->chain->rev_seeds),
        kv_size(data->chain->rev_seeds),
        FALSE,
        REV,
        &hit_list,
        &hit_seed_list);
    s_fill_init_hit_info(data->query_oid, data->query_size,
        data->subject_oid, data->subject_size,
        kv_data(hit_list), kv_size(hit_list));
    for (size_t i = 0; i < kv_size(hit_list); ++i) {
        HbnInitHit* hit = &kv_A(hit_list, i);
        hit->chain_seed_offset = kv_size(data->hit_seed_list);
        kv_push(HbnInitHit, data->hit_list, *hit);
        kv_push_v(ChainSeed, data->hit_seed_list, hit->chain_seed_array, hit->chain_seed_count);
    }
    kv_destroy(hit_list);
    kv_destroy(hit_seed_list);

    for (size_t i = 0; i < kv_size(data->hit_list); ++i) {
        kv_A(data->hit_list, i).chain_seed_array = kv_data(data->hit_seed_list) 
                                                   + 
                                                   kv_A(data->hit_list, i).chain_seed_offset;
    }
    ks_introsort_init_hit_score_gt(kv_size(data->hit_list), kv_data(data->hit_list));
}

void
InitHitFindData_SetupMapAlignInfoFromHit(InitHitFindData* data,
    HbnInitHit* hit,
    vec_chain_seed* seed_list)
{
    hbn_assert(hit->qdir == FWD);
    hbn_assert(hit->sdir == FWD || hit->sdir == REV);
    const u8* q_in_seed = data->fwd_query;
    hbn_assert(q_in_seed != NULL);
    const u8* s_in_seed = (hit->sdir == FWD) ? data->fwd_subject : data->rev_subject;
    hbn_assert(s_in_seed != NULL);

    kv_clear(*seed_list);
    int i = 0;
    while (i < hit->chain_seed_count) {
        ChainSeed* ci = (ChainSeed*)hit->chain_seed_array + i;
        ChainSeed* cj = NULL;
        int j = i + 1;
        while (j < hit->chain_seed_count) {
            cj = (ChainSeed*)hit->chain_seed_array + j;
            if (cj->qoff - ci->qoff != cj->soff - ci->soff) break;
            ++j;            
        }
        cj = (ChainSeed*)hit->chain_seed_array + (j - 1);
        int ml = cj->qoff + cj->length - ci->qoff;
        hbn_assert(ml > 0);       
        s_extract_mem(s_in_seed, q_in_seed,
            ci->soff, ci->qoff, ml, data->word_size, seed_list);
        i = j;
    }

    validate_mem(HBN_LOG_ARGS_DEFAULT,
        s_in_seed,
        q_in_seed,
        kv_data(*seed_list),
        kv_size(*seed_list));
}

void
InitHitFindData_SetupCnsAlignInfoFromHit(InitHitFindData* data,
    HbnInitHit* hit,
    vec_chain_seed* seed_list)
{
    hbn_assert(hit->qdir == FWD);
    hbn_assert(hit->sdir == FWD || hit->sdir == REV);
    const u8* q_in_seed = data->fwd_query;
    hbn_assert(q_in_seed != NULL);
    const u8* s_in_seed = (hit->sdir == FWD) ? data->fwd_subject : data->rev_subject;
    hbn_assert(s_in_seed != NULL);

    kv_clear(*seed_list);
    int i = 0;    
    while (i < hit->chain_seed_count){
        ChainSeed* ci =  (ChainSeed*)hit->chain_seed_array + i;
        ChainSeed* cj = NULL;
        int j = i + 1;
        while (j < hit->chain_seed_count) {
            cj = (ChainSeed*)hit->chain_seed_array + j;
            if (cj->qoff - ci->qoff != cj->soff - ci->soff) break;
            ++j;
        }
        cj = (ChainSeed*)hit->chain_seed_array + (j - 1);
        int ml = cj->qoff + cj->length - ci->qoff;
        s_extract_mem(q_in_seed, s_in_seed,
            ci->qoff, ci->soff, ml, data->word_size, seed_list);
        i = j;
    }
    if (kv_empty(*seed_list)) return;

    if (hit->sdir == REV) {
        ChainSeed* csa = kv_data(*seed_list);
        int csc = kv_size(*seed_list);
        for (int i = 0; i < csc; ++i) {
            hbn_assert(csa[i].qoff + csa[i].length <= data->query_size);
            int qoff = data->query_size - csa[i].qoff - csa[i].length;
            hbn_assert(csa[i].soff + csa[i].length <= data->subject_size);
            int soff = data->subject_size - csa[i].soff - csa[i].length;
            csa[i].qoff = qoff;
            csa[i].soff = soff;
        }

        int s = 0, e = csc - 1;
        while (s < e) {
            ChainSeed tmp = csa[s];
            csa[s] = csa[e];
            csa[e] = tmp;
            ++s;
            --e;
        }
    }

    validate_mem(HBN_LOG_ARGS_DEFAULT,
        (hit->sdir == FWD) ? data->fwd_query : data->rev_query,
        data->fwd_subject,
        kv_data(*seed_list),
        kv_size(*seed_list));
}