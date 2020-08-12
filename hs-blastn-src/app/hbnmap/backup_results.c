#include "backup_results.h"

#include "traceback_stage.h"
#include "../../ncbi_blast/setup/hsp2string.h"

static int
extract_query_id(BlastHitList* hit_list)
{
    int query_id = -1;
    for (int i = 0; i < hit_list->hsplist_count; ++i) {
        BlastHSPList* hsp_list = hit_list->hsplist_array[i];
        if (!hsp_list) continue;
        for (int j = 0; j < hsp_list->hspcnt; ++j) {
            BlastHSP* hsp = hsp_list->hsp_array[j];
            if (!hsp) continue;
            int oid = hsp->hbn_query.oid;
            query_id = oid;
            break;
        }
    }
    return query_id;
}

static void
add_one_hsp(const BlastHSP* hsp, kstring_t* out)
{
    const char* in = (const char*)(hsp);
    int in_len = sizeof(BlastHSP);
    kputsn(in, in_len, out);

    int num_op = hsp->gap_info ? hsp->gap_info->size : 0;
    in = (const char*)(&num_op);
    in_len = sizeof(int);
    kputsn(in, in_len, out);
    if (!num_op) return;

    in = (const char*)(hsp->gap_info->op_type);
    in_len = sizeof(EGapAlignOpType) * num_op;
    kputsn(in, in_len, out);

    in = (const char*)(hsp->gap_info->num);
    in_len = sizeof(int) * num_op;
    kputsn(in, in_len, out);
}

static BlastHSP*
read_one_hsp(const char** in_)
{
    const char* in = *in_;
    BlastHSP* hsp = Blast_HSPNew();
    memcpy(hsp, in, sizeof(BlastHSP));
    in += sizeof(BlastHSP);
    int num_op = 0;
    memcpy(&num_op, in, sizeof(int));
    in += sizeof(int);
    if (num_op) {
        hsp->gap_info = GapEditScriptNew(num_op);
        int len = sizeof(EGapAlignOpType) * num_op;
        memcpy(hsp->gap_info->op_type, in, len);
        in += len;
        len = sizeof(int) * num_op;
        memcpy(hsp->gap_info->num, in, len);
        in += len;
    }
    *in_ = in;
    return hsp;
}

static BlastHSPList*
read_one_hsp_list(const char** in_)
{
    const char* in = *in_;
    int num_hsp = 0;
    memcpy(&num_hsp, in, sizeof(int));
    in += sizeof(int);
    BlastHSPList* hsp_list = (BlastHSPList*)calloc(1, sizeof(BlastHSPList));
    hsp_list->hsp_array = (BlastHSP**)calloc(num_hsp, sizeof(BlastHSP*));
    hsp_list->hspcnt = num_hsp;
    hsp_list->hsp_max = num_hsp;
    for (int i = 0; i < num_hsp; ++i) hsp_list->hsp_array[i] = read_one_hsp(&in);
    hsp_list->oid = hsp_list->query_index = hsp_list->hsp_array[0]->hbn_query.oid;
    double best_evalue = hsp_list->hsp_array[0]->evalue;
    double best_raw_score = hsp_list->hsp_array[0]->score;
    for (int i = 0; i < num_hsp; ++i) {
        best_evalue = hbn_min(best_evalue, hsp_list->hsp_array[i]->evalue);
        best_raw_score = hbn_max(best_raw_score, hsp_list->hsp_array[i]->score);
    }
    hsp_list->best_evalue = best_evalue;
    hsp_list->hbn_best_raw_score = best_raw_score;
    *in_ = in;
    return hsp_list;
}

static void
add_one_hsp_list(BlastHSPList* hsp_list, kstring_t* out)
{
    if (hsp_list->hspcnt == 0) return;
    const char* in = (const char*)(&hsp_list->hspcnt);
    int len = sizeof(int);
    kputsn(in, len, out);
    for (int i = 0; i < hsp_list->hspcnt; ++i) {
        add_one_hsp(hsp_list->hsp_array[i], out);
    }
}

static void
add_one_hit_list(BlastHitList* hit_list, kstring_t* out)
{
    if (hit_list->hsplist_count == 0) return;
    const char* in = (const char*)(&hit_list->hsplist_count);
    int len = sizeof(int);
    kputsn(in, len, out);
    for (int i = 0; i < hit_list->hsplist_count; ++i) {
        add_one_hsp_list(hit_list->hsplist_array[i], out);
    }
}

static void
read_one_hit_list(BlastHitList* hit_list, const char** in_)
{
    const char* in = *in_;
    int hsplist_count = 0;
    memcpy(&hsplist_count, in, sizeof(int));
    in += sizeof(int);
    hit_list->hsplist_array = (BlastHSPList**)calloc(hsplist_count, sizeof(BlastHSPList*));
    for (int i = 0; i < hsplist_count; ++i) hit_list->hsplist_array[i] = read_one_hsp_list(&in);

    hit_list->hsplist_count = hit_list->hsplist_max = hsplist_count;
    double worst_evalue = 0.0;
    int low_score = INT32_MAX;
    int num_hits = 0;
    for (int i = 0; i < hsplist_count; ++i) {
        BlastHSPList* hsp_list = hit_list->hsplist_array[i];
        for (int j = 0; j < hsp_list->hspcnt; ++j) {
            BlastHSP* hsp = hsp_list->hsp_array[j];
            ++num_hits;
            worst_evalue = hbn_max(worst_evalue, hsp->evalue);
            low_score = hbn_min(low_score, hsp->score);
        }
    }
    hit_list->num_hits = num_hits;
    hit_list->low_score = low_score;
    hit_list->worst_evalue = worst_evalue;
    *in_ = in;
}

BOOL 
read_one_hsp_result_set(FILE* in, HbnHSPResults* results)
{
    int len = 0;
    int n = fread(&len, sizeof(int), 1, in);
    if (n == 0) return FALSE;
    HbnHSPResultsClear(results, 0);
    char* ins = (char*)calloc(len, sizeof(char));
    n = fread(ins, 1, len, in);
    hbn_assert(n == len);
    char* ins_end = ins + len;
    const char* p = ins;
    int num_queries = 0;
    while (p < ins_end) {
        read_one_hit_list(&results->hitlist_array[num_queries], &p);
        ++num_queries;
    }
    results->num_queries = num_queries;
    free(ins);
    return TRUE;
}

void
add_one_hsp_result_set(HbnHSPResults* results, FILE* out, pthread_mutex_t* out_lock)
{
    kstring_t* os = &results->output_buf;
    ks_clear(*os);
    if (results->num_queries == 0) return;
    for (int i = 0; i < results->num_queries; ++i) add_one_hit_list(&results->hitlist_array[i], os);
    int len = ks_size(*os);
    if (out_lock) pthread_mutex_lock(out_lock);
    hbn_fwrite(&len, sizeof(int), 1, out);
    hbn_fwrite(ks_s(*os), 1, ks_size(*os), out);
    if (out_lock) pthread_mutex_unlock(out_lock);
    ks_clear(*os);
}

void
dump_m4_hits(const text_t* query_vol,
    const text_t* subject_vol,
    HbnHSPResults* results,
    const HbnProgramOptions* opts)
{
    ks_dinit(line);
    ks_clear(results->output_buf);
    for (int i = 0; i < results->num_queries; ++i) {
        BlastHitList* hit_list = results->hitlist_array + i;
        if (hit_list->hsplist_array == 0) continue;
        for (int j = 0; j < hit_list->hsplist_count; ++j) {
            BlastHSPList* hsp_list = hit_list->hsplist_array[j];
            for (int k = 0; k < hsp_list->hspcnt; ++k) {
                //HBN_LOG("*** hspcnt = %d", hsp_list->hspcnt);
                BlastHSP* hsp = hsp_list->hsp_array[k];
                hbn_assert(hsp->hbn_subject.strand == FWD);
                if (hsp->hbn_query.strand == REV) {
                    int offset = hsp->hbn_query.seq_size - hsp->hbn_query.end;
                    int end = hsp->hbn_query.seq_size - hsp->hbn_query.offset;
                    hsp->hbn_query.offset = offset;
                    hsp->hbn_query.end = end;
                }
                if (hsp->hbn_subject.strand == REV) {
                    int offset = hsp->hbn_subject.seq_size - hsp->hbn_subject.end;
                    int end = hsp->hbn_subject.seq_size - hsp->hbn_subject.offset;
                    hsp->hbn_subject.offset = offset;
                    hsp->hbn_subject.end = end;
                }
                const char* qname = seqdb_seq_name(query_vol, hsp->hbn_query.oid);
                const char* sname = seqdb_seq_name(subject_vol, hsp->hbn_subject.oid);
                hsp->hbn_query.oid = query_vol->dbinfo.seq_start_id + hsp->hbn_query.oid;
                hsp->hbn_subject.oid = subject_vol->dbinfo.seq_start_id + hsp->hbn_subject.oid;

                print_one_sam_result(hsp,
                    &results->aligned_strings,
                    qname,
                    sname,
                    opts->dump_md,
                    opts->rg_sample,
                    &results->output_buf);
            }
        }
    }
    ks_destroy(line);
}

static void
recover_aligned_strings(const CSeqDB* queries,
    const CSeqDB* db,
    HbnHSPResults* results)
{
    kv_dinit(vec_u8, fwd_query);
    kv_dinit(vec_u8, rev_query);
    for (int i = 0; i < results->num_queries; ++i) {
        BlastHitList* hit_list = results->hitlist_array + i;
        int query_id = extract_query_id(hit_list);
        if (query_id == -1) continue;
        seqdb_extract_sequence(queries, query_id, FWD, &fwd_query);
        seqdb_recover_sequence_ambig_res(queries, query_id, FWD, kv_data(fwd_query));
        seqdb_extract_sequence(queries, query_id, REV, &rev_query);
        seqdb_recover_sequence_ambig_res(queries, query_id, REV, kv_data(rev_query));
        for (int j = 0; j < hit_list->hsplist_count; ++j) {
            BlastHSPList* hsp_list = hit_list->hsplist_array[j];
            if (hsp_list->hspcnt == 0) continue;
            const int subject_id = hsp_list->hsp_array[0]->hbn_subject.oid;
            const u8* subject = db->unpacked_seq + seqdb_seq_offset(db, subject_id);
            for (int k = 0; k < hsp_list->hspcnt; ++k) {
                BlastHSP* hsp = hsp_list->hsp_array[k];
                const u8* query = (hsp->hbn_query.strand == FWD) ? kv_data(fwd_query) : kv_data(rev_query);
                add_align_string(hsp, query, subject, &results->aligned_strings);
            }
        }
    }
    kv_destroy(fwd_query);
    kv_destroy(rev_query);
}

static void
purge_null_hsplist(HbnHSPResults* results)
{
    for (int i = 0; i < results->num_queries; ++i) {
        BlastHitList* hit_list = results->hitlist_array + i;
        int n = 0;
        //HBN_LOG("hsplist: %d", hit_list->hsplist_count);
        for (int j = 0; j < hit_list->hsplist_count; ++j) {
            BlastHSPList* hsp_list = hit_list->hsplist_array[j];
            if (hsp_list) {
                Blast_HSPListPurgeNullHSPs(hsp_list);
                if (hsp_list->hspcnt == 0) {
                    free(hsp_list);
                    hsp_list = NULL;
                }
            }
            if (hsp_list) hit_list->hsplist_array[n++] = hsp_list;
        }
        //HBN_LOG("n = %d", n);
        hit_list->hsplist_count = n;
    }
}

void
dump_one_result_set(const CSeqDB* queries,
    const CSeqDB* db,
    HbnHSPResults* results, 
    const HbnProgramOptions* opts,
    FILE* out, 
    FILE* backup_out, 
    pthread_mutex_t* out_lock)
{
    purge_null_hsplist(results);

    if (opts->outfmt == eSAM) {
        recover_aligned_strings(queries, db, results);
        dump_m4_hits(queries, db, results, opts);
    } else if (opts->outfmt == eTabular || opts->outfmt == eTabularWithComments) {
        print_tabular_reports(results, opts->subject, db, queries, opts->outfmt);
    }
    if (out_lock) pthread_mutex_lock(out_lock);
    hbn_fwrite(ks_s(results->output_buf), 1, ks_size(results->output_buf), out);
    if (out_lock) pthread_mutex_unlock(out_lock);

    if (backup_out) add_one_hsp_result_set(results, backup_out, out_lock);
}

void
recover_qi_vs_sj_results(const CSeqDB* queries, const CSeqDB* db, const HbnProgramOptions* opts, FILE* in, FILE* out)
{
    HbnHSPResults* results = HbnHSPResultsNew(HBN_QUERY_CHUNK_SIZE);
    while (read_one_hsp_result_set(in, results)) {
        dump_one_result_set(queries, db, results, opts, out, NULL, NULL);
    }
    results = HbnHSPResultsFree(results);
}