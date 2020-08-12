#include "primer_map_one_volume.h"

#include "blast_gapalign.h"
#include "primer_map_hit_finder.h"
#include "../hbnmap/backup_results.h"
#include "../hbnmap/search_setup.h"
#include "../hbnmap/tabular_format.h"
#include "../hbnmap/traceback_stage.h"
#include "../../algo/hbn_traceback_aux.h"
#include "../../ncbi_blast/setup/blast_encoding.h"
#include "../../ncbi_blast/setup/blast_hits.h"
#include "../../ncbi_blast/setup/hsp2string.h"
#include "../../ncbi_blast/setup/blast_query_info.h"
#include "../../ncbi_blast/setup/blast_sequence_blk.h"

#include <pthread.h>
#include <stdio.h>

static FILE* g_out = NULL;
static pthread_mutex_t g_out_lock;
static int g_qid = 0;
static pthread_mutex_t g_qid_lock;
static const HbnProgramOptions* g_opts = NULL;
static const HbnOptionsHandle* g_opts_handle = NULL;
static CSeqDB* g_primer_volume = NULL;
static CSeqDB* g_query_volume = NULL;

static const int kPrimerBatchSize = 100;
static const int kQueryBatchSize = 100;

static void
s_init_global_data(CSeqDB* p_vol, CSeqDB* q_vol, 
    const HbnProgramOptions* opts,
    const HbnOptionsHandle* opts_handle,
    FILE* out)
{
    g_out = out;
    pthread_mutex_init(&g_out_lock, NULL);
    g_qid = 0;
    pthread_mutex_init(&g_qid_lock, NULL);
    g_opts = opts;
    g_opts_handle = opts_handle;
    g_primer_volume = p_vol;
    g_query_volume = q_vol;
}

static BOOL
s_extract_query_batch(PrimerMapHitFindData* hit_data,
    BLAST_SequenceBlk* query_blk,
    BlastQueryInfo* query_info)
{
    //HBN_LOG("load query");
    hbn_assert(g_query_volume != NULL);
    if (!extract_sequence_block_from_packed_seqdb(g_query_volume,
            &g_qid,
            &g_qid_lock,
            kQueryBatchSize,
            TRUE,
            TRUE,
            TRUE,
            query_blk,
            query_info)) return FALSE;
    PrimerMapHitFindData_CleanQueryData(hit_data);
    for (int i = 0; i < query_info->num_queries; ++i) {
        int ctx_id = i * 2;
        int oid = query_info->contexts[ctx_id].query_index;
        const u8* query = query_blk->sequence + query_info->contexts[ctx_id].query_offset;
        int query_length = query_info->contexts[ctx_id].query_length;
        PrimerMapHitFindData_AddOneQuery(hit_data,
            oid,
            query,
            query_length,
            seqdb_seq_name(g_query_volume, oid));
    }

    return TRUE;
}

static BOOL
s_extract_primer_batch(PrimerMapHitFindData* hit_data, 
    int* pid,
    pthread_mutex_t* pid_lock,
    BLAST_SequenceBlk* primer_blk,
    BlastQueryInfo* primer_info)
{
    //HBN_LOG("load primer");
    if (!extract_sequence_block_from_unpacked_seqdb(g_primer_volume,
            pid,
            pid_lock,
            kPrimerBatchSize,
            TRUE,
            TRUE,
            TRUE,
            primer_blk,
            primer_info)) return FALSE;

    PrimerMapHitFindData_CleanPrimerData(hit_data);
    for (int i = primer_info->first_context; i <= primer_info->last_context; ++i) {
        int ctx_id = i;
        int oid = primer_info->contexts[ctx_id].query_index;
        const u8* primer = primer_blk->sequence + primer_info->contexts[ctx_id].query_offset;
        int primer_length = primer_info->contexts[ctx_id].query_length;
        PrimerMapHitFindData_AddOnePrimer(hit_data,
            oid,
            primer,
            primer_length,
            seqdb_seq_name(g_primer_volume, oid));
    }
    return TRUE;    
}

static void
s_extract_mem(const u8* primer,
    const u8* subject,
    const int poff,
    const int soff,
    const int length,
    int* p_start,
    int* s_start,
    int* max_l)
{
    int i = 0;
    int pi = poff, si = soff;
    while (i < length) {
        while (i < length) {
            if (primer[pi] == subject[si]) break;
            ++pi;
            ++si;
            ++i;
        }
        if (i >= length) break;

        hbn_assert(primer[pi] == subject[si]);
        int j = i;
        while (j < length) {
            if (primer[pi] != subject[si]) break;
            ++pi;
            ++si;
            ++j;
        }

        int len = j - i;
        if (len >= (*max_l)) {
            *p_start = pi - len / 2;
            *s_start = si - len / 2;
            *max_l = len;
        }
        i = j;
    }
}

static void
s_find_gapped_start(const u8* primer,
    const u8* subject,
    HbnInitHit* hit, int* p_start, int* s_start)
{
    int max_l = 0;
    ChainSeed* csa = (ChainSeed*)hit->chain_seed_array;
    int csc = hit->chain_seed_count;
    int i = 0;
    while (i < csc) {
        int j = i + 1;
        while (j < csc) {
            if (csa[j].qoff - csa[i].qoff != csa[j].soff - csa[i].soff) break;
            ++j;
        }
        int len = csa[j-1].soff + csa[j-1].length - csa[i].soff;
        s_extract_mem(primer, subject, csa[i].qoff, csa[i].soff, len, p_start, s_start, &max_l);
        i = j;
    }
}

static BOOL
s_seed_is_contained_in_hsp_list(BlastHSP* hsp_list, int hsp_cnt, int qdir, int qoff, int poff)
{
    const int E = 2;
    for (int i = 0; i < hsp_cnt; ++i) {
        const BlastHSP* hsp = hsp_list + i;
        hbn_assert(hsp->hbn_subject.strand == FWD);
        if (hsp->hbn_query.strand != qdir) continue;
        int r = (qoff + E >= hsp->hbn_query.offset)
                &&
                (qoff <= hsp->hbn_query.end + E)
                &&
                (poff + E >= hsp->hbn_subject.offset)
                &&
                (poff <= hsp->hbn_subject.end + E);
        if (r) return TRUE;
    }
    return FALSE;
}

static void 
set_blasthsp(BlastGapAlignStruct* gap_align,
    int qid,
    const int fwd_query_context,
    int qdir,
    int qsize,
    int sid,
    int subject_offset,
    int ssize,
    const char* qaln,
    const char* saln,
    int aln_size,
    const HbnProgramOptions* opts,
    int ddf_score,
    int chain_score,
    BlastHSP* hsp)
{
    int score = 0;
    int reward = abs(opts->reward);
    int penalty = -abs(opts->penalty);
    int go = -abs(opts->gap_open);
    int ge = -abs(opts->gap_extend);
    int i = 0;
    int num_op = 0;
    int num_ident = 0;
    while (i < aln_size) {
        char qc = qaln[i];
        char sc = saln[i];
        if (qc != GAP_CHAR && sc != GAP_CHAR) {
            int j = i;
            while (j < aln_size) {
                qc = qaln[j];
                sc = saln[j];
                if (qc == GAP_CHAR || sc == GAP_CHAR) break;
                score += (qc == sc) ? reward : penalty;
                num_ident += (qc == sc);
                ++j;
            }
            ++num_op;
            i = j;
            continue;
        }
        if (qc == GAP_CHAR) {
            int l = 0;
            int j = i;
            while (j < aln_size) {
                qc = qaln[j];
                sc = saln[j];
                if (qc != GAP_CHAR) break;
                ++l;
                ++j;
            }
            score += (go + ge * l);
            ++num_op;
            i = j;
            continue;
        }
        hbn_assert(sc == GAP_CHAR);
        {
            int l = 0;
            int j = i;
            while (j < aln_size) {
                qc = qaln[j];
                sc = saln[j];
                if (sc != GAP_CHAR) break;
                ++l;
                ++j;
            }
            score += (go + ge * l);
            ++num_op;
            i = j;
            continue;
        }
    }

    hsp->score = score;
    hsp->num_ident = num_ident;
    hsp->bit_score = 0.0;
    hsp->evalue = 0.0;
    hsp->query.frame = qdir;
    hsp->query.offset = gap_align->query_start;
    hsp->query.end = gap_align->query_stop;
    hsp->query.gapped_start = gap_align->query_start;
    hsp->subject.frame = 0;
    hsp->subject.offset = subject_offset + gap_align->subject_start;
    hsp->subject.end = subject_offset + gap_align->subject_stop;
    hsp->subject.gapped_start = hsp->subject.offset;
    hsp->context = fwd_query_context + qdir;
    hsp->num_positives = num_ident;

    hsp->gap_info = GapEditScriptNew(num_op);
    int op_idx = 0;
    int gap_opens = 0;
    int gaps = 0;
    i = 0;
    while (i < aln_size) {
        char qc = qaln[i];
        char sc = saln[i];
        if (qc != GAP_CHAR && sc != GAP_CHAR) {
            int j = i;
            EGapAlignOpType type = eGapAlignSub;
            int l = 0;
            while (j < aln_size) {
                qc = qaln[j];
                sc = saln[j];
                if (qc == GAP_CHAR || sc == GAP_CHAR) break;
                ++l;
                ++j;
            }
            hsp->gap_info->op_type[op_idx] = type;
            hsp->gap_info->num[op_idx] = l;
            ++op_idx;
            i = j;
            continue;
        }
        if (qc == GAP_CHAR) {
            int l = 0;
            int j = i;
            EGapAlignOpType type = eGapAlignDel;
            while (j < aln_size) {
                qc = qaln[j];
                sc = saln[j];
                if (qc != GAP_CHAR) break;
                ++l;
                ++j;
            }
            hsp->gap_info->op_type[op_idx] = type;
            hsp->gap_info->num[op_idx] = l;
            ++gap_opens;
            gaps += l;
            ++op_idx;
            i = j;
            continue;
        }
        hbn_assert(sc == GAP_CHAR);
        {
            int l = 0;
            int j = i;
            EGapAlignOpType type = eGapAlignIns;
            while (j < aln_size) {
                qc = qaln[j];
                sc = saln[j];
                if (sc != GAP_CHAR) break;
                ++l;
                ++j;
            }
            hsp->gap_info->op_type[op_idx] = type;
            hsp->gap_info->num[op_idx] = l;
            ++gap_opens;
            gaps += l;
            ++op_idx;
            i = j;
            continue;
        }
    }    

    hsp->hbn_query.oid = qid;
    hsp->hbn_query.strand = qdir;
    hsp->hbn_query.offset = hsp->query.offset;
    hsp->hbn_query.end = hsp->query.end;
    hsp->hbn_query.seq_size = qsize;

    hsp->hbn_subject.oid = sid;
    hsp->hbn_subject.strand = FWD;
    hsp->hbn_subject.offset = hsp->subject.offset;
    hsp->hbn_subject.end = hsp->subject.end;
    hsp->hbn_subject.seq_size = ssize;

    hsp->hsp_info.perc_identity = 100.0 * num_ident / aln_size;
    hsp->hsp_info.ddf_score = ddf_score;
    hsp->hsp_info.chain_score = chain_score;
    hsp->hsp_info.num_ident = num_ident;
    hsp->hsp_info.num_positives = num_ident;
    hsp->hsp_info.align_len = aln_size;
    hsp->hsp_info.gap_opens = gap_opens;
    hsp->hsp_info.gaps = gaps;

    if (op_idx != num_op) {
        HBN_LOG("op_idx = %d, num_op = %d, align_size = %d", op_idx, num_op, aln_size);
        dump_align_string(qaln, saln, aln_size, stderr);
        dump_blasthsp(fprintf, stderr, *hsp);
    }
    hbn_assert(op_idx == num_op);
}

void
s_extract_align_string_from_ges(int qoff, int qend, int soff, int send, 
    const u8* query, const u8* subject, GapEditScript* gap_info,
    kstring_t* qaln,
    kstring_t* saln)
{
    const char BLASTNA_TO_IUPACNA[BLASTNA_SIZE] = {
    'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 
    'W', 'S', 'B', 'D', 'H', 'V', 'N', '-'
    };
    ks_clear(*qaln);
    ks_clear(*saln);
    int qi = qoff;
    const u8* q = query + qi;
    for (int i = 0; i < gap_info->size; ++i) {
        EGapAlignOpType type = gap_info->op_type[i];
        int num = gap_info->num[i];
        if (type == eGapAlignSub || type == eGapAlignIns) {
            for (int p = 0; p < num; ++p, ++q, ++qi) {
                u8 c = *q;
                int dc = BLASTNA_TO_IUPACNA[c];
                kputc(dc, qaln);
            }
        } else if (type == eGapAlignDel) {
            for (int p = 0; p < num; ++p) {
                kputc(GAP_CHAR, qaln);
            }
        }
    }
    hbn_assert(qi == qend);

    int si = soff;
    const u8* s = subject + si;
    for (int i = 0; i < gap_info->size; ++i) {
        EGapAlignOpType type = gap_info->op_type[i];
        int num = gap_info->num[i];
        if (type == eGapAlignSub || type == eGapAlignDel) {
            for (int p = 0; p < num; ++p, ++s, ++si) {
                u8 c = *s;
                int dc = BLASTNA_TO_IUPACNA[c];
                kputc(dc, saln);
            }
        } else if (type == eGapAlignIns) {
            for (int p = 0; p < num; ++p) {
                kputc(GAP_CHAR, saln);
            }
        }
    }
    hbn_assert(si == send);
    hbn_assert(ks_size(*qaln) == ks_size(*saln));
}

static void
qx_map_one_query_subject_hit_list(PrimerMapHitFindData* hit_finder,
    BlastGapAlignStruct* gap_align,
    BlastScoringParameters* score_params,
    BLAST_SequenceBlk* query_blk,
    BlastQueryInfo* query_info,
    BLAST_SequenceBlk* primer_blk,
    BlastQueryInfo* primer_info,
    HbnInitHit* hit_array,
    int hit_count,
    BlastHSPList* hsp_list)
{
    BlastHSP hsp_array[g_opts->max_hsps_per_subject];
    int hspcnt = 0;
    ks_introsort_init_hit_score_gt(hit_count, hit_array);
    ks_dinit(qaln);
    ks_dinit(saln);
    int primer_index = kv_A(hit_finder->primer_context_list, hit_array[0].sid).seq_index;
    int query_index = kv_A(hit_finder->query_context_list, hit_array[0].qid).seq_index;
    const char* query_name = kv_A(hit_finder->query_context_list, hit_array[0].qid).name;
    //if (strcmp(query_name, "03381ecc-1c0f-4128-bfec-0a30e740e631")) return;
    //HBN_LOG("number of hits: %d", hit_count);
    for (int i = 0; i < hit_count && i < g_opts->max_hsps_per_subject + 5; ++i) {
        HbnInitHit* hit = hit_array + i;
        hbn_assert(hit->qdir == FWD);
        hbn_assert(hit->sdir == FWD);
        const u8* primer = PrimerMapHitFindData_ExtractPrimer(hit_finder, hit->sid);
        const int primer_size = hit->ssize;
        int primer_dir = FWD;
        const u8* query = PrimerMapHitFindData_ExtractQuery(hit_finder, hit->qid);
        const int query_size = hit->qsize;
        int query_dir = FWD;
        int q_start, p_start;
        s_find_gapped_start(query, primer, hit, &q_start, &p_start);
        //HBN_LOG("q_start = %d, p_start = %d", q_start, p_start);

        if (hit->sid & 1) {
            query_dir = REV;
            query = query_blk->sequence
                    +
                    query_info->contexts[hit->qid*2+1].query_offset;
            q_start = query_size - 1 - q_start;
            hbn_assert(q_start >= 0);
            primer_dir = FWD;
            primer = primer_blk->sequence
                     +
                     primer_info->contexts[hit->sid-1].query_offset;
            p_start = primer_size - 1 - p_start;
            hbn_assert(p_start >= 0);
        }

        if (s_seed_is_contained_in_hsp_list(hsp_array, hspcnt, query_dir, q_start, p_start)) continue;
#if 0
        Int2 status = BLAST_GappedAlignmentWithTraceback(eBlastTypeBlastn,
                        query,
                        primer,
                        gap_align,
                        score_params,
                        q_start,
                        p_start,
                        query_size,
                        primer_size,
                        NULL);
#else
        Int2 status = BLAST_GreedyGappedAlignment(query, primer,
                        query_size, primer_size, gap_align,
                        score_params, 
                        q_start, p_start,
                        FALSE, TRUE, NULL);
#endif 
        if (status != 0) continue;
        s_extract_align_string_from_ges(gap_align->query_start,
            gap_align->query_stop,
            gap_align->subject_start,
            gap_align->subject_stop,
            query,
            primer,
            gap_align->edit_script,
            &qaln,
            &saln);
        //dump_align_string(ks_s(qaln), ks_s(saln), ks_size(qaln), stderr);
        BlastHSP* hsp = hsp_array + hspcnt;
        memset(hsp, 0, sizeof(BlastHSP));
        ++hspcnt;
        set_blasthsp(gap_align,
            query_index,
            hit->qid * 2,
            query_dir,
            query_size,
            primer_index,
            0,
            primer_size,
            ks_s(qaln),
            ks_s(saln),
            ks_size(qaln),
            g_opts,
            hit->score,
            hit->score,
            hsp);
        //dump_blasthsp(fprintf, stderr, *hsp);
        if (hspcnt == g_opts->max_hsps_per_subject) break;
    }
    //HBN_LOG("hspcnt = %d", hspcnt);
    ks_destroy(qaln);
    ks_destroy(saln);
    if (hspcnt == 0) return;

    ks_introsort_blasthsp_score_gt(hspcnt, hsp_array);
    hsp_list->oid = primer_index;
    hsp_list->hbn_best_raw_score = hsp_array[0].score;
    hsp_list->hspcnt = hspcnt;
    hsp_list->hsp_max = hspcnt;
    hsp_list->hsp_array = (BlastHSP**)calloc(hspcnt, sizeof(BlastHSP*));
    for (int i = 0; i < hspcnt; ++i) {
        BlastHSP* hsp = Blast_HSPNew();
        memcpy(hsp, hsp_array + i, sizeof(BlastHSP));
        hsp_list->hsp_array[i] = hsp;
    }
    //exit(0);
}

static void
qx_map_one_query_hit_list(PrimerMapHitFindData* hit_finder,
    BlastGapAlignStruct* gap_align,
    BlastScoringParameters* score_params,
    BLAST_SequenceBlk* query_blk,
    BlastQueryInfo* query_info,
    BLAST_SequenceBlk* primer_blk,
    BlastQueryInfo* primer_info,    
    HbnInitHit* hit_array,
    int hit_count,
    BlastHitList* hit_list)
{
    BlastHSPList hsplist_array[kPrimerBatchSize];
    memset(hsplist_array, 0, sizeof(BlastHSPList) * kPrimerBatchSize);
    int hsplist_count = 0;
    int i = 0;
    while (i < hit_count) {
        int max_sid = hit_array[i].sid / 2;
        max_sid = max_sid * 2 + 1;
        int j = i + 1;
        while (j < hit_count && hit_array[j].sid <= max_sid) ++j;
        qx_map_one_query_subject_hit_list(hit_finder,
            gap_align,
            score_params,
            query_blk,
            query_info,
            primer_blk,
            primer_info,
            hit_array + i,
            j - i,
            &hsplist_array[hsplist_count]);
        if (hsplist_array[hsplist_count].hspcnt > 0) ++hsplist_count;
        hbn_assert(hsplist_count <= kPrimerBatchSize);
        i = j;
    }

    if (!hsplist_count) return;
    hit_list->hsplist_array = (BlastHSPList**)calloc(hsplist_count, sizeof(BlastHSPList*));
    hit_list->hsplist_count = hsplist_count;
    hit_list->hsplist_max = hsplist_count;
    for (i = 0; i < hsplist_count; ++i) {
        BlastHSPList* hsp_list = (BlastHSPList*)calloc(1, sizeof(BlastHSPList));
        memcpy(hsp_list, hsplist_array + i, sizeof(BlastHSPList));
        hit_list->hsplist_array[i] = hsp_list;
    }
}

static void
qx_map_one_batch(PrimerMapHitFindData* hit_finder,
    BLAST_SequenceBlk* primer_blk,
    BlastQueryInfo* primer_info,
    BLAST_SequenceBlk* query_blk,
    BlastQueryInfo* query_info,
    HbnHSPResults* results)
{
    BlastScoreBlk* sbp = NULL;
    BlastEffectiveLengthsParameters* eff_len_params = NULL;
    BlastScoringParameters* score_params = NULL;
    BlastExtensionParameters* ext_params = NULL;
    BlastHitSavingParameters* hit_params = NULL;
    BlastInitialWordParameters* word_params = NULL;
    const HbnOptionsHandle* opts_handle = g_opts_handle;
    sbp = CSetupFactory__CreateScoreBlock(opts_handle, query_blk, query_info);
    BlastScoreBlkCheck(sbp);
    BLAST_GapAlignSetUp(eBlastTypeBlastn,
        g_primer_volume,
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
    BlastGapAlignStruct* gap_align = NULL;
    BLAST_GapAlignStructNew(score_params,
        ext_params,
        MAX_DBSEQ_LEN,
        sbp,
        &gap_align);
    HbnHSPResultsClear(results, query_info->num_queries);

    PrimerMapHitFindData_FindHits(hit_finder);
    HbnInitHit* hit_array = kv_data(hit_finder->hit_list);
    int hit_count = kv_size(hit_finder->hit_list);
    int i = 0;
    for (i = 0; i < hit_count - 1; ++i) hbn_assert(hit_array[i].qid <= hit_array[i+1].qid);
    i = 0;
    //HBN_LOG("number of hits: %d", hit_count);
    while (i < hit_count) {
        const int qid = hit_array[i].qid;
        int j = i + 1;
        while (j < hit_count && hit_array[j].qid == qid) ++j;
        hbn_assert(qid < results->num_queries);
        //HBN_LOG("qid = %d", qid);
        BlastHitList* hit_list = results->hitlist_array + qid;
        qx_map_one_query_hit_list(hit_finder,
            gap_align,
            score_params,
            query_blk,
            query_info,
            primer_blk,
            primer_info,
            hit_array + i,
            j - i,
            hit_list);
        i = j;

        //HBN_LOG("number of hsplist: %d", hit_list->hsplist_count);
        for (j = 0; j < hit_list->hsplist_count; ++j) {
            compute_traceback_from_hsplist(eBlastTypeBlastn, 
                hit_list->hsplist_array[j], 
                query_blk, 
                query_info, 
                g_primer_volume, 
                sbp, 
                score_params, 
                ext_params->options, 
                hit_params,
                &results->aligned_strings);
            //HBN_LOG("qid = %d, sid = %d, hspcnt = %d", 
            //    qid, hit_list->hsplist_array[j]->oid, hit_list->hsplist_array[j]->hspcnt);
        }
        //HBN_LOG("number of hsplist: %d", hit_list->hsplist_count);
    }
    dump_one_result_set(g_query_volume, g_primer_volume, results, g_opts, g_out, NULL, &g_out_lock);

    sbp = BlastScoreBlkFree(sbp);
    word_params = BlastInitialWordParametersFree(word_params);
    hit_params = BlastHitSavingParametersFree(hit_params);
    ext_params = BlastExtensionParametersFree(ext_params);
    score_params = BlastScoringParametersFree(score_params);
    eff_len_params = BlastEffectiveLengthsParametersFree(eff_len_params);
    gap_align = BLAST_GapAlignStructFree(gap_align);
}

static void*
qx_map_thread(void* params)
{
    PrimerMapHitFindData* hit_finder = PrimerMapHitFindDataNew(g_opts->memsc_kmer_size,
                                            g_opts->memsc_kmer_window,
                                            g_opts->memsc_score);
    HbnHSPResults* results = HbnHSPResultsNew(kQueryBatchSize);
    BLAST_SequenceBlk* primer_blk = BLAST_SequenceBlkNew();
    BlastQueryInfo* primer_info = BlastQueryInfoNew(2 * kPrimerBatchSize);
    BLAST_SequenceBlk* query_blk = BLAST_SequenceBlkNew();
    BlastQueryInfo* query_info = BlastQueryInfoNew(2 * kQueryBatchSize);
    pthread_mutex_t pid_lock = PTHREAD_MUTEX_INITIALIZER;

    while (s_extract_query_batch(hit_finder, query_blk, query_info)) {
        PrimerMapHitFindData_BuildQueryWordList(hit_finder);
        int pid = 0;
        while (s_extract_primer_batch(hit_finder, &pid, &pid_lock, primer_blk, primer_info)) {
            PrimerMapHitFindData_BuildPrimerWordList(hit_finder);
            qx_map_one_batch(hit_finder, primer_blk, primer_info, query_blk, query_info, results);
        }
        //break;
    }

    hit_finder = PrimerMapHitFindDataFree(hit_finder);
    results = HbnHSPResultsFree(results);
    primer_blk = BLAST_SequenceBlkFree(primer_blk);
    primer_info = BlastQueryInfoFree(primer_info);
    query_blk = BLAST_SequenceBlkFree(query_blk);
    query_info = BlastQueryInfoFree(query_info);
    return NULL;
}

void
primer_map_one_volume(CSeqDB* qvol, 
    CSeqDB* svol, 
    const HbnProgramOptions* opts,
    const HbnOptionsHandle* opts_handle,
    FILE* out)
{
    s_init_global_data(svol, qvol, opts, opts_handle, out);
    pthread_t jobs[opts->num_threads];
    for (int i = 0; i < opts->num_threads; ++i) {
        pthread_create(jobs + i, NULL, qx_map_thread, NULL);
    }
    for (int i = 0; i < opts->num_threads; ++i) {
        pthread_join(jobs[i], NULL);
    }
}