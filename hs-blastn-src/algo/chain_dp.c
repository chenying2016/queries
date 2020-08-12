#include "chain_dp.h"

#include "../corelib/ksort.h"

void
validate_mem(HBN_LOG_PARAMS_GENERIC,
    const u8* read, 
    const u8* subject,
    const ChainSeed* cdpsa,
    const int cdpsc)
{
    //return;
    //HBN_LOG("validating mems from [%s, %s, %d]", HBN_LOG_ARGS_GENERIC);
    for (int i = 0; i < cdpsc; ++i) {
        ChainSeed s = cdpsa[i];
        //fprintf(stderr, "\tvalidating %d\t%d\t%d\t%d\t%d\n", i, s.qoff, s.soff, s.length, s.sdir);
        int qi = s.qoff;
        int si = s.soff;
        for (int k = 0; k < s.length; ++k, ++qi, ++si) {
            hbn_assert(read[qi] == subject[si], "[%s, %s, %d] at (%d, %d, %d): qi = %d, si = %d", 
            HBN_LOG_ARGS_GENERIC, s.qoff, s.soff, s.length, qi, si);
        }
    }
}

ChainWorkData*
ChainWorkDataNew(int min_seed_cnt, int min_can_score)
{
    ChainWorkData* data = (ChainWorkData*)calloc(1, sizeof(ChainWorkData));
    kv_init(data->f);
    kv_init(data->p);
    kv_init(data->t);
    kv_init(data->v);
    kv_init(data->u);
    kv_init(data->seeds);
    kv_init(data->fwd_seeds);
    kv_init(data->rev_seeds);

    data->max_dist_qry = 5000;
    data->max_dist_ref = 5000;
    data->max_band_width = 1500;

#if 1
    data->max_join_long = 20000;
    data->max_join_short = 2000;
    data->min_join_flank_sc = 1000;
    data->min_join_flank_ratio = 0.5;
#else 
    data->max_join_long = 1000;
    data->max_join_short = 200;
    data->min_join_flank_sc = 200;
    data->min_join_flank_ratio = 0.5;
#endif

    data->max_skip = 25;
    data->min_cnt = min_seed_cnt;
    data->min_score = min_can_score;
    return data;
}

ChainWorkData*
ChainWorkDataFree(ChainWorkData* data)
{
    kv_destroy(data->f);
    kv_destroy(data->p);
    kv_destroy(data->t);
    kv_destroy(data->v);
    kv_destroy(data->u);
    kv_destroy(data->seeds);
    kv_destroy(data->fwd_seeds);
    kv_destroy(data->rev_seeds);
    free(data); 
    return NULL;   
}

static void
ChainWorkDataSetup(ChainWorkData* data, int n)
{
    kv_resize(int, data->f, n);
    kv_resize(int, data->p, n);
    kv_resize(int, data->t, n);
    kv_resize(int, data->v, n);
    kv_zero(int, data->t);
    kv_zero(int, data->f);
    kv_fill(data->p, -1);
    kv_resize(IntPair, data->u, n);
}

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

static void
scoring_chain_seeds(ChainWorkData* data, 
    const ChainSeed* seeds,
    const int n,
    const BOOL is_maximal_exact_match)
{
    const int max_dist_ref = data->max_dist_ref;
    const int max_dist_qry = data->max_dist_qry;
    const int band_width = data->max_band_width;
    const int max_skip = data->max_skip;
    int sum_cov = 0;
    for (int i = 0; i < n; ++i) sum_cov += seeds[i].length;
    const int avg_cov = sum_cov / n;
    int st = 0;
    ChainWorkDataSetup(data, n);
    int* f = kv_data(data->f);
    int* p = kv_data(data->p);
    int* t = kv_data(data->t);
    int* v = kv_data(data->v);

    // fill the score and backtrack arrays
    for (int i = 0; i < n; ++i) {
        idx ri = seeds[i].soff;
        idx qi = seeds[i].qoff;
        int max_j = -1;
        int cov = seeds[i].length;
        int max_f = cov, n_skip = 0, min_d;
        while (st < i && ri > seeds[st].soff + seeds[st].length + max_dist_ref) ++st;
        for (int j = i - 1; j >= st; --j) {
            //HBN_LOG("comparing");
            //fprintf(stderr, "[%d, %d, %d, %d] v.s. [%d, %d, %d, %d]\n",
            //    i, seeds[i].qoff, seeds[i].soff, seeds[i].length, j, seeds[j].qoff, seeds[j].soff, seeds[j].length);
            if (is_maximal_exact_match) {
                if (seeds[j].qoff + seeds[j].length > qi || seeds[j].soff + seeds[j].length > ri) continue;
            } else {
                if (seeds[j].qoff > qi || seeds[j].soff > ri) continue;
            }
            idx dr = ri - seeds[j].soff;
            idx dq = qi - seeds[j].qoff;
            int dd, sc, log_dd;
            if (dr == 0 || dq <= 0) continue;
            if (dq > max_dist_qry || dr > max_dist_ref) continue;
            dd = (dr > dq) ? (dr - dq) : (dq - dr);
            if (dd > band_width) continue;
            min_d = hbn_min(dq, dr);
            sc = (min_d > cov) ? cov : hbn_min(dq, dr);
            log_dd = dd ? ilog2_32(dd) : 0;
            if (!is_maximal_exact_match) sc -= (int)(dd * .01 * avg_cov) + (log_dd>>1);
            else sc -= (log_dd>>1);
            sc += f[j];
            if (sc > max_f) {
                max_f = sc;
                max_j = j;
                if (n_skip) --n_skip;
            } else if (t[j] == i) {
                if (++n_skip > max_skip) { break; }
            }
            if (p[j] >= 0) t[p[j]] = i;
        }
        f[i] = max_f;
        p[i] = max_j;
        //HBN_LOG("f[%d] = %d, p[%d] = %d", i, max_f, i, p[i]);
        // v[i] keeps the peak score up to i;
        // f[i] is the score ending at i, not always the peak score
        v[i] = (max_j >= 0 && v[max_j] > max_f) ? v[max_j] : max_f;
    }
}

int
find_best_kmer_match(ChainWorkData* data,
    int* best_kmer_match_index,
    int* best_kmer_match_score)
{
    const int n = kv_size(data->seeds);
    if (n == 0) return 0;
    const int min_cnt = data->min_cnt;
    const int min_score = data->min_score;
    scoring_chain_seeds(data, kv_data(data->seeds), kv_size(data->seeds), FALSE);
    int* f = kv_data(data->f);
    int* p = kv_data(data->p);
    int* t = kv_data(data->t);
    int* v = kv_data(data->v);
    IntPair* u = kv_data(data->u);

    // find the ending position of chains
    memset(t, 0, sizeof(int) * n);
    for (int i = 0; i < n; ++i) {
        if (p[i] >= 0) t[p[i]] = 1;
    }
    int n_u;
    for (int i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_score) ++n_u;
    }
    if (n_u == 0) return 0;
    for (int i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_score) {
            int j = i;
            while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maxizes f
            if (j < 0) {
                j = i;
            }
            u[n_u].first = f[j];
            u[n_u].second = j;
            ++n_u;
        }
    }
    ks_introsort_int_pair(n_u, u);
    // reverse u, such that highest scoring chains appear first
    for (int i = 0; i < n_u>>1; ++i) {
        IntPair tmp = u[i];
        u[i] = u[n_u - i - 1];
        u[n_u - i - 1] = tmp;
    }

    // backtrack
    memset(t, 0, sizeof(int) * n);
    int n_v, k;
    int find_best_seed = 0;
    for (int i = n_v = k = 0; i < n_u; ++i) { // start from the highest score
        int n_v0 = n_v;
        int j = u[i].second;
        do {
            v[n_v++] = j;
            t[j] = 1;
            j = p[j];
        } while (j >= 0 && t[j] == 0);
        if (j < 0) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first;
                u[k].second = n_v - n_v0;
                ++k;
                find_best_seed = 1;
            }
        } else if (u[i].first - f[j] >= min_score) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first - f[j];
                u[k].second = n_v - n_v0;
                ++k;
                find_best_seed = 1;
            }
        }

        if (find_best_seed) {
            if (best_kmer_match_index) *best_kmer_match_index = v[n_v0];
            if (best_kmer_match_score) *best_kmer_match_score = u[k-1].first;
            return 1;
        } else {
            return 0;
        }
    }    
    return 0;
}

static BOOL
init_hit_is_contained(HbnInitHit* ha, int hc, HbnInitHit* hn)
{
    const int E = 100;
    for (int i = 0; i < hc; ++i) {
        int r = (hn->qbeg + E >= ha[i].qbeg)
                &&
                (hn->qend <= ha[i].qend + E)
                &&
                (hn->sbeg + E >= ha[i].sbeg)
                &&
                (hn->send <= ha[i].send + E);
        if (r) return TRUE;
    }
    return FALSE;
}

int chaining_find_candidates(ChainWorkData* data,
        ChainSeed* chain_seed_array,
        int chain_seed_count,
        const int is_maximal_exact_match,
        const int subject_strand,
        vec_init_hit* init_hit_list,
        vec_chain_seed* chain_seed_list)
{
    const ChainSeed* seeds = chain_seed_array;
    const int n = chain_seed_count;
    if (n == 0) return 0;
    const int min_cnt = data->min_cnt;
    const int min_score = data->min_score;
    //HBN_LOG("socring mems");
    scoring_chain_seeds(data, chain_seed_array, chain_seed_count, is_maximal_exact_match);
    //HBN_LOG("done");
    int* f = kv_data(data->f);
    int* p = kv_data(data->p);
    int* t = kv_data(data->t);
    int* v = kv_data(data->v);
    IntPair* u = kv_data(data->u);
    HbnInitHit hit = { 0, FWD, 0, FWD, 0 };
    hit.qdir = FWD;
    hit.sdir = subject_strand;

    // find the ending position of chains
    memset(t, 0, sizeof(int) * n);
    for (int i = 0; i < n; ++i) {
        if (p[i] >= 0) t[p[i]] = 1;
    }
    int n_u;
    for (int i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_score) ++n_u;
    }
    if (n_u == 0) return 0;
    for (int i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_score) {
            int j = i;
            while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maxizes f
            if (j < 0) {
                j = i;
            }
            u[n_u].first = f[j];
            u[n_u].second = j;
            ++n_u;
        }
    }
    ks_introsort_int_pair(n_u, u);
    // reverse u, such that highest scoring chains appear first
    for (int i = 0; i < n_u>>1; ++i) {
        IntPair tmp = u[i];
        u[i] = u[n_u - i - 1];
        u[n_u - i - 1] = tmp;
    }

    // backtrack
    memset(t, 0, sizeof(int) * n);
    int n_v, k;
    int find_can = 0;
    int num_added_hit = 0;
    int first_hit_idx = kv_size(*init_hit_list);
    for (int i = n_v = k = 0; i < n_u && num_added_hit < 40; ++i) { // start from the highest score
        int n_v0 = n_v, k0 = k;
        int j = u[i].second;
        if (t[j]) continue;
        find_can = 0;
        do {
            v[n_v++] = j;
            t[j] = 1;
            j = p[j];
        } while (j >= 0 && t[j] == 0);
        if (j < 0) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first;
                u[k].second = n_v - n_v0;
                ++k;
                find_can = 1;
            }
        } else if (u[i].first - f[j] >= min_score) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first - f[j];
                u[k].second = n_v - n_v0;
                ++k;
                find_can = 1;
            }
        }

        if (find_can) {
            hit.score = u[k-1].first;
            hit.chain_seed_offset = kv_size(*chain_seed_list);
            hit.chain_seed_count = n_v - n_v0;
            int max_size = 0;
            hit.qbeg = seeds[v[n_v-1]].qoff;
            hit.sbeg = seeds[v[n_v-1]].soff;
            hit.qend = seeds[v[n_v0]].qoff + seeds[v[n_v0]].length;
            hit.send = seeds[v[n_v0]].soff + seeds[v[n_v0]].length;
            if (!init_hit_is_contained(kv_data(*init_hit_list) + first_hit_idx, num_added_hit, &hit)) {
                for (int x = n_v; x > n_v0; --x) {
                    int y = v[x-1];
                    kv_push(ChainSeed, *chain_seed_list, seeds[y]);
                    if (seeds[y].length > max_size) {
                        max_size = seeds[y].length;
                        hit.qoff = seeds[y].qoff + max_size / 2;
                        hit.soff = seeds[y].soff + max_size / 2;
                    }
                }
		//HBN_LOG("find init hit, [%d, %d] x [%d, %d], score = %d", hit.qbeg, hit.qend, hit.sbeg, hit.send, hit.score);
                //dump_init_hit(fprintf, stderr, hit);
                kv_push(HbnInitHit, *init_hit_list, hit);
                ++num_added_hit;
            }
        }

        if (k0 == k) n_v = n_v0;
    }    
    join_adjacent_chains(data, init_hit_list, chain_seed_list);
    for (size_t p = 0; p < kv_size(*init_hit_list); ++p) {
        HbnInitHit* hit = &kv_A(*init_hit_list, p);
        hit->chain_seed_array = kv_data(*chain_seed_list) + hit->chain_seed_offset;
        hit->chain_seed_offset = INT32_MAX;
    }
    return 0;
}

static void
set_init_hit_pos_from_chain_seed_list(HbnInitHit* hit, ChainSeed* chain_seed_array)
{
    ChainSeed cs = chain_seed_array[hit->chain_seed_offset];
    hit->qbeg = cs.qoff;
    hit->sbeg = cs.soff;
    cs = chain_seed_array[hit->chain_seed_offset + hit->chain_seed_count-1];
    hit->qend = cs.qoff + cs.length;
    hit->send = cs.soff + cs.length;
}

#define init_hit_soff_lt(a, b) (((a).sbeg < (b).sbeg) || ((a).sbeg == (b).sbeg && (a).qbeg < (b).qbeg))
KSORT_INIT(init_hit_soff_lt, HbnInitHit, init_hit_soff_lt);

#define intpair_gt(a, b) (((a).first > (b).first) || ((a).first == (b).first && (a).second < (b).second))
KSORT_INIT(intpair_gt, IntPair, intpair_gt);

static BOOL
two_chains_are_adjacent(ChainWorkData* data, HbnInitHit* left, HbnInitHit* right)
{
    if (left->qend > right->qbeg || left->send > right->sbeg) return FALSE;
    int gap_qry = right->qbeg - left->qend;
    int gap_ref = right->sbeg - left->send;
    int max_gap = hbn_max(gap_qry, gap_ref);
    int min_gap = hbn_min(gap_qry, gap_ref);
    if (max_gap > data->max_join_long || min_gap > data->max_join_short) return FALSE;

    int left_qlen = left->qend - left->qbeg;
    int left_slen = left->send - left->sbeg;
    int right_qlen = right->qend - right->qbeg;
    int right_slen = right->send - left->sbeg;
    const int kMinMemLen = 1000;
    const int kMinMemScore = 500;
    int r = (left_qlen < kMinMemLen) || (right_qlen < kMinMemLen) || (left_slen < kMinMemLen) || (right_slen < kMinMemLen);
    if (r) return FALSE;
    r = (left->score < kMinMemScore) || (right->score < kMinMemScore);
    if (r) return FALSE;
    return TRUE;

#if 0
    int sc_thres = (int)(1.0 * data->min_join_flank_sc / data->max_join_long * max_gap + .499);
    if (left->score < sc_thres || right->score < sc_thres) return FALSE;
    int min_flank_len = (int)(1.0 * max_gap * data->min_join_flank_ratio);
    if (left->send - left->sbeg < min_flank_len || left->qend - left->qbeg < min_flank_len) return FALSE;
    if (right->send - right->sbeg < min_flank_len || right->qend - right->qbeg < min_flank_len) return FALSE;
    return TRUE;
#endif
}

void
join_adjacent_chains(ChainWorkData* data, vec_init_hit* init_hit_list, vec_chain_seed* chain_seed_list)
{
    if (kv_size(*init_hit_list) < 2) return;

    const int kMaxExamineInitHit = 20;
    HbnInitHit* hit_array = kv_data(*init_hit_list);
    int hit_count = kv_size(*init_hit_list);
    for (int i = 0; i < hit_count; ++i) {
        set_init_hit_pos_from_chain_seed_list(hit_array + i, kv_data(*chain_seed_list));
    }
    ks_introsort_init_hit_soff_lt(hit_count, hit_array);
    IntPair* hit_score_and_idx_array = (IntPair*)calloc(hit_count, sizeof(IntPair));
    for (int i = 0; i < hit_count; ++i) {
        IntPair ip = { hit_array[i].score, i };
        hit_score_and_idx_array[i] = ip;
    }
    ks_introsort_intpair_gt(hit_count, hit_score_and_idx_array);

    HbnInitHit new_hit_array[kMaxExamineInitHit];
    int new_hit_count = 0;
    kv_dinit(vec_chain_seed, new_chain_seed_list);
    kv_dinit(vec_init_hit, hit_chain);
    for (int i = 0; i < hit_count && i < kMaxExamineInitHit; ++i) {
        int hit_idx = hit_score_and_idx_array[i].second;
        HbnInitHit* hit = hit_array + hit_idx;
        if (hit->qdir == 2) continue;
        HbnInitHit new_hit = *hit;
        hit->qdir = 2;
        int k = hit_idx;
        kv_clear(hit_chain);
        kv_push(HbnInitHit, hit_chain, new_hit);
        while (k) {
            --k;
            hit = hit_array + k;
            if (hit->qdir == 2) continue;
            if (!two_chains_are_adjacent(data, hit, &new_hit)) continue;
            kv_push(HbnInitHit, hit_chain, *hit);
            new_hit.qbeg = hit->qbeg;
            new_hit.sbeg = hit->sbeg;
            new_hit.score += hit->score;
            hit->qdir = 2;
        }
        k = hit_idx;
        while (k < hit_count - 1) {
            ++k;
            hit = hit_array + k;
            if (hit->qdir == 2) continue;
            if (!two_chains_are_adjacent(data, &new_hit, hit)) continue;
            kv_push(HbnInitHit, hit_chain, *hit);
            new_hit.qend = hit->qend;
            new_hit.send = hit->send;
            new_hit.score += hit->score;
            hit->qdir = 2;
        }

        ks_introsort_init_hit_soff_lt(kv_size(hit_chain), kv_data(hit_chain));
        new_hit.chain_seed_offset = kv_size(new_chain_seed_list);
        new_hit.chain_seed_count = 0;
        for (size_t p = 0; p < kv_size(hit_chain); ++p) {
            hit = &kv_A(hit_chain, p);
            ChainSeed* csa = kv_data(*chain_seed_list) + hit->chain_seed_offset;
            int csc = hit->chain_seed_count;
            new_hit.chain_seed_count += csc;
            kv_push_v(ChainSeed, new_chain_seed_list, csa, csc);
        }
        new_hit_array[new_hit_count++] = new_hit;
#if 0
        if (kv_size(hit_chain) > 2) {
            HBN_LOG("merge %zu hits", kv_size(hit_chain));
            for (size_t p = 0; p < kv_size(hit_chain); ++p) {
                HbnInitHit* hit1 = &kv_A(hit_chain, p);   
                fprintf(stderr, "\t[%d, %d] x [%d, %d], score = %d, offset = %d, cnt = %d\n",
                    hit1->qbeg, hit1->qend, hit1->sbeg, hit1->send, hit1->score, 
                    hit1->chain_seed_offset, hit1->chain_seed_count);
            }
        }
#endif 
    }

    ks_introsort_init_hit_score_gt(new_hit_count, new_hit_array);
    kv_clear(*init_hit_list);
    kv_push_v(HbnInitHit, *init_hit_list, new_hit_array, new_hit_count);
    kv_clear(*chain_seed_list);
    kv_push_v(ChainSeed, *chain_seed_list, kv_data(new_chain_seed_list), kv_size(new_chain_seed_list));
    kv_destroy(new_chain_seed_list);
    kv_destroy(hit_chain);
    free(hit_score_and_idx_array);
}
