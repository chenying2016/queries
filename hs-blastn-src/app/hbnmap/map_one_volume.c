#include "map_one_volume.h"

#include "../../algo/chain_dp.h"
#include "../../algo/hbn_traceback.h"
#include "../../algo/init_hit_finder.h"
#include "../../algo/word_finder.h"

#include "../../corelib/gapped_candidate.h"
#include "../../corelib/seqdb.h"

#include "../../ncbi_blast/setup/blast_parameters.h"
#include "../../ncbi_blast/setup/blast_query_info.h"
#include "../../ncbi_blast/setup/blast_sequence_blk.h"
#include "../../ncbi_blast/setup/hsp2string.h"

#include "search_setup.h"
#include "subseq_hit.h"
#include "find_seeding_subseqs.h"
#include "hbn_find_subseq_hit.h"
#include "traceback_stage.h"
#include "hbn_results.h"
#include "tabular_format.h"
#include "backup_results.h"

#include <limits.h>
#include <pthread.h>

static int g_thread_id = -1;
static pthread_mutex_t g_thread_id_lock;
static int g_query_index = -1;
static pthread_mutex_t g_query_index_lock;

static void
init_global_data()
{
    g_thread_id = 0;
    pthread_mutex_init(&g_thread_id_lock, NULL);
    g_query_index = 0;
    pthread_mutex_init(&g_query_index_lock, NULL);
}

static int
get_next_query_chunk(
    const text_t* queries, 
    int* next_query_id,
    pthread_mutex_t* query_id_lock,
    BLAST_SequenceBlk* query_blk, 
    BlastQueryInfo* query_info)
{
    int from = 0, to = 0;
    pthread_mutex_lock(query_id_lock);
    from = *next_query_id;
    *next_query_id += HBN_QUERY_CHUNK_SIZE;
    pthread_mutex_unlock(query_id_lock);
    if (from >= queries->dbinfo.num_seqs) return 0;
    to = hbn_min(from + HBN_QUERY_CHUNK_SIZE, queries->dbinfo.num_seqs);
    int num_queries = to - from;

    int length = 0;
    for (int i = from; i < to; ++i) length += seqdb_seq_size(queries, i);
    length *= 2;

    BlastContextInfo ctx_info;
    int ctx_idx = 0;
    int seq_idx = 0;
    int max_length = 0;
    int min_length = I32_MAX;
    query_blk->sequence = (Uint1*)realloc(query_blk->sequence, length);
    query_blk->sequence_nomask = (Uint1*)realloc(query_blk->sequence_nomask, length);
    query_blk->length = length;
    kv_dinit(vec_u8, seq);

    for (int i = from; i < to; ++i) {
        ctx_info.query_length = seqdb_seq_size(queries, i);
        ctx_info.eff_searchsp = 0;
        ctx_info.length_adjustment = 0;
        ctx_info.query_index = i;
        ctx_info.frame = 0;
        ctx_info.is_valid = TRUE;
        ctx_info.segment_flags = 0;

        max_length = hbn_max(max_length, ctx_info.query_length);
        min_length = hbn_min(min_length, ctx_info.query_length);

        ctx_info.query_offset = seq_idx;
        query_info->contexts[ctx_idx++] = ctx_info;
        seqdb_extract_sequence(queries, i, FWD, &seq);
        hbn_assert(ctx_info.query_length == kv_size(seq));
        memcpy(query_blk->sequence + seq_idx, kv_data(seq), kv_size(seq));
        seqdb_recover_sequence_ambig_res(queries, i, FWD, kv_data(seq));
        hbn_assert(ctx_info.query_length == kv_size(seq));
        memcpy(query_blk->sequence_nomask + seq_idx, kv_data(seq), kv_size(seq));
        seq_idx += ctx_info.query_length;

        ctx_info.query_offset = seq_idx;
        query_info->contexts[ctx_idx++] = ctx_info;
        seqdb_extract_sequence(queries, i, REV, &seq);
        hbn_assert(ctx_info.query_length == kv_size(seq));
        memcpy(query_blk->sequence + seq_idx, kv_data(seq), kv_size(seq));
        seqdb_recover_sequence_ambig_res(queries, i, REV, kv_data(seq));
        hbn_assert(ctx_info.query_length == kv_size(seq));
        memcpy(query_blk->sequence_nomask + seq_idx, kv_data(seq), kv_size(seq));
        seq_idx += ctx_info.query_length;
    }
    hbn_assert(seq_idx == length);

    query_info->first_context = 0;
    query_info->last_context = ctx_idx - 1;
    query_info->num_queries = num_queries;
    query_info->max_length = max_length;
    query_info->min_length = min_length;  

    kv_destroy(seq);
    return num_queries;  
}

static void
reverse_seeding_subseqs(vec_int_pair* seeding_subseq_list, const int query_length)
{
    if (kv_empty(*seeding_subseq_list)) return;

    IntPair* seeding_subseq_array = kv_data(*seeding_subseq_list);
    int seeding_subseq_count = kv_size(*seeding_subseq_list);
    IntPair ip;
    for (int i = 0; i < seeding_subseq_count; ++i) {
        ip.first = query_length - seeding_subseq_array[i].second;
        ip.second = query_length - seeding_subseq_array[i].first;
        seeding_subseq_array[i] = ip;
    }
    int left = 0, right = seeding_subseq_count - 1;
    while (left < right) {
        ip = seeding_subseq_array[left];
        seeding_subseq_array[left] = seeding_subseq_array[right];
        seeding_subseq_array[right] = ip;
        ++left;
        --right;
    }
}

void
align_one_query_block(CSeqDB* query_vol,
    CSeqDB* subject_vol,
    WordFindData* word_data,
    HbnSubseqHitExtnData* extn_data,
    BLAST_SequenceBlk* query_blk,
    BlastQueryInfo* query_info,
    const HbnProgramOptions* opts,
    const HbnOptionsHandle* opts_handle,
    HbnHSPResults* results,
    FILE* out,
    FILE* backup_out,
    pthread_mutex_t* out_lock)
{
    BlastScoreBlk* sbp = NULL;
    BlastEffectiveLengthsParameters* eff_len_params = NULL;
    BlastScoringParameters* score_params = NULL;
    BlastExtensionParameters* ext_params = NULL;
    BlastHitSavingParameters* hit_params = NULL;
    BlastInitialWordParameters* word_params = NULL;
    sbp = CSetupFactory__CreateScoreBlock(opts_handle, query_blk, query_info);
    BlastScoreBlkCheck(sbp);
    BLAST_GapAlignSetUp(eBlastTypeBlastn,
        subject_vol,
        opts_handle->m_ScoringOpts,
        opts_handle->m_EffLenOpts,
        opts_handle->m_ExtnOpts,
        opts_handle->m_HitSaveOpts,
        opts_handle->m_InitWordOpts,
        query_info,
        sbp,
        &score_params,
        &ext_params,
        &hit_params,
        &eff_len_params,
        &word_params);
    HbnHSPResultsClear(results, query_info->num_queries);

    vec_int_pair* seeding_subseq_list = &word_data->seeding_subseqs;
    kv_dinit(vec_subseq_hit, fwd_subseq_hit_list);
    kv_dinit(vec_subseq_hit, rev_subseq_hit_list);
    kv_dinit(vec_subseq_hit, subseq_hit_list);
    kv_dinit(vec_u8, subject_v);
    
    for (int i = 0; i < query_info->num_queries; ++i) {
        int ctx_id = i * 2;
        BlastContextInfo ctx_info = query_info->contexts[ctx_id];
        int query_id = ctx_info.query_index;
        int query_strand = FWD;
        int query_length = ctx_info.query_length;
        const char* query_name = seqdb_seq_name(query_vol, query_id);
        const u8* fwd_query = query_blk->sequence + ctx_info.query_offset;
        kv_clear(word_data->init_hit_list);
        kv_clear(*seeding_subseq_list);
        if (opts->strand == FWD || opts->strand == F_R) {
            find_seeding_subseqs(fwd_query, query_length, opts->kmer_size, seeding_subseq_list);
            ddfs_find_candidates(word_data, fwd_query, query_id, query_vol->dbinfo.seq_start_id, query_strand, query_length);
        }

        ctx_id++;
        ctx_info = query_info->contexts[ctx_id];
        query_strand = REV;
        const u8* rev_query = query_blk->sequence + ctx_info.query_offset;
        if (opts->strand == REV || opts->strand == F_R) {
            if (kv_empty(*seeding_subseq_list)) {
                find_seeding_subseqs(rev_query, query_length, opts->kmer_size, seeding_subseq_list);
            } else {
                reverse_seeding_subseqs(seeding_subseq_list, query_length);
            }
            ddfs_find_candidates(word_data, rev_query, query_id, query_vol->dbinfo.seq_start_id, query_strand, query_length);
        }

        HbnInitHit* init_hit_array = kv_data(word_data->init_hit_list);
        int init_hit_count = kv_size(word_data->init_hit_list);
        find_candidate_subject_subseqs(init_hit_array,
            init_hit_count,
            query_id,
            query_length,
            opts,
            subject_vol,
            &fwd_subseq_hit_list,
            &rev_subseq_hit_list,
            &subseq_hit_list);

        hbn_extend_query_subseq_hit_list(kv_data(subseq_hit_list),
            kv_size(subseq_hit_list),
            query_name,
            fwd_query,
            rev_query,
            query_id,
            ctx_id - 1,
            query_length,
            &subject_v,
            subject_vol,
            opts,
            extn_data,
            results->hitlist_array + i,
            results);
        
        BlastHitList* hit_list = results->hitlist_array + i;
        for (int j = 0; j < hit_list->hsplist_count; ++j) {
            compute_traceback_from_hsplist(eBlastTypeBlastn, 
                hit_list->hsplist_array[j], 
                query_blk, 
                query_info, 
                subject_vol, 
                sbp, 
                score_params, 
                ext_params->options, 
                hit_params,
                &results->aligned_strings);
        }
    }

    dump_one_result_set(query_vol, subject_vol, results, opts, out, backup_out, out_lock);

    kv_destroy(fwd_subseq_hit_list);
    kv_destroy(rev_subseq_hit_list);
    kv_destroy(subseq_hit_list);
    kv_destroy(subject_v);
    sbp = BlastScoreBlkFree(sbp);
    word_params = BlastInitialWordParametersFree(word_params);
    hit_params = BlastHitSavingParametersFree(hit_params);
    ext_params = BlastExtensionParametersFree(ext_params);
    score_params = BlastScoringParametersFree(score_params);
    eff_len_params = BlastEffectiveLengthsParametersFree(eff_len_params);

    int gid = query_info->contexts[0].query_index;
    if (gid && (gid % 1000 == 0)) HBN_LOG("%8d queries processed", gid);
}

static void*
hbn_align_worker(void* params)
{
    hbn_task_struct* ht_struct = (hbn_task_struct*)(params);
    int thread_id = -1;
    pthread_mutex_lock(&g_thread_id_lock);
    thread_id = g_thread_id;
    ++g_thread_id;
    pthread_mutex_unlock(&g_thread_id_lock);
    hbn_assert(thread_id >= 0 && thread_id < ht_struct->opts->num_threads);
    CSeqDB* query_vol = ht_struct->query_vol;
    CSeqDB* subject_vol = ht_struct->subject_vol;
    WordFindData* word_data = ht_struct->word_data_array[thread_id];
    HbnSubseqHitExtnData* extn_data = ht_struct->hit_extn_data_array[thread_id];
    const HbnProgramOptions* opts = ht_struct->opts;
    const HbnOptionsHandle* opts_handle = ht_struct->opts_handle;
    HbnHSPResults* results = ht_struct->results_array[thread_id];
    BLAST_SequenceBlk* query_blk = BLAST_SequenceBlkNew();
    BlastQueryInfo* query_info = BlastQueryInfoNew(HBN_QUERY_CHUNK_SIZE * 2);

    while (get_next_query_chunk(query_vol,
                &g_query_index,
                &g_query_index_lock,
                query_blk,
                query_info)) {
        align_one_query_block(query_vol,
            subject_vol,
            word_data,
            extn_data,
            query_blk,
            query_info,
            opts,
            opts_handle,
            results,
            ht_struct->out,
            ht_struct->qi_vs_sj_out,
            &ht_struct->out_lock);
    }

    BLAST_SequenceBlkFree(query_blk);
    BlastQueryInfoFree(query_info);
    return NULL;
}

void
hbn_align_one_volume(hbn_task_struct* ht_struct)
{
    init_global_data();
    const int num_threads = ht_struct->opts->num_threads;
    pthread_t job_ids[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(job_ids + i, NULL, hbn_align_worker, ht_struct);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(job_ids[i], NULL);
    }
}