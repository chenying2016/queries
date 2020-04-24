#include "hbn_extend_subseq_hit.h"

#include "../../algo/hbn_traceback_aux.h"
#include "../../ncbi_blast/setup/hsp2string.h"

HbnSubseqHitExtnData*
HbnSubseqHitExtnDataNew(const HbnProgramOptions* opts)
{
    HbnSubseqHitExtnData* data = (HbnSubseqHitExtnData*)calloc(1, sizeof(HbnSubseqHitExtnData));
        const int memsc_kmer_size = opts->memsc_kmer_size;
        const int memsc_mem_size = opts->memsc_mem_size;
        hbn_assert(memsc_kmer_size <= memsc_mem_size);
        const int memsc_window_size = opts->memsc_kmer_window;
        const int memsc_score = opts->memsc_score;
        data->mem_data = MaximalExactMatchWorkDataNew(memsc_kmer_size, 
                            memsc_window_size, 
                            memsc_mem_size, 
                            memsc_score);
        data->traceback_data = HbnTracebackDataNew();
    kv_init(data->chain_seed_list);
    kv_init(data->fwd_sbjct_subseq_list);
    kv_init(data->rev_sbjct_subseq_list);
    kv_init(data->sbjct_subseq_list);
    return data;
}

HbnSubseqHitExtnData*
HbnSubseqHitExtnDataFree(HbnSubseqHitExtnData* data)
{
    if (data->mem_data) MaximalExactMatchWorkDataFree(data->mem_data);
    if (data->traceback_data) HbnTracebackDataFree(data->traceback_data);
    kv_destroy(data->chain_seed_list);
    kv_destroy(data->fwd_sbjct_subseq_list);
    kv_destroy(data->rev_sbjct_subseq_list);
    kv_destroy(data->sbjct_subseq_list);
    free(data);
    return NULL;
}

static BOOL
subseq_is_contained(const BlastHSP* hsp_array,
    const int hsp_count,
    HbnSubseqHit* hit,
    const int qsize)
{
    const int E = 100;
    int qoff = hit->qoff;
    int soff = hit->soff + hit->sfrom;
    for (int i = 0; i < hsp_count; ++i) {
        const BlastHSP* hsp = hsp_array + i;
        if (hsp->hbn_query.strand != hit->qdir) continue;
        int r = (qoff + E >= hsp->hbn_query.offset)
                &&
                (qoff <= hsp->hbn_query.end + E)
                &&
                (soff + E >= hsp->hbn_subject.offset)
                &&
                (soff <= hsp->hbn_subject.end + E);
        if (r) return TRUE;
    }

    int hit_qoff = hit->qoff;
    if (hit->qdir == REV) hit_qoff = qsize - 1 - hit->qoff;
    int hit_soff = hit->soff + hit->sfrom;
    for (int i = -1; i < hsp_count; ++i) {
        const BlastHSP* hsp = (hsp_array + i);
        if (hsp->hsp_info.ddf_score == -1) continue;
        int hsp_qbeg = hsp->hbn_query.offset;
        int hsp_qend = hsp->hbn_query.end;
        if (hsp->hbn_query.strand == REV) {
            hsp_qbeg = qsize - hsp->hbn_query.end;
            hsp_qend = qsize - hsp->hbn_query.offset;
        }
        int hsp_sbeg = hsp->hbn_subject.offset;
        int hsp_send = hsp->hbn_subject.end;
        int r = (hit->qoff >= hsp_qbeg && hit->qoff <= hsp_qend)
                ||
                (hit_soff >= hsp_sbeg && hit_soff <= hsp_send);
        if (!r) continue;
        if (hit->score < hsp->hsp_info.ddf_score * 0.4) return TRUE;
    }
    return FALSE;
}

static BOOL
chain_seed_list_is_contained(const BlastHSP* hsp_array,
    const int hsp_count,
    const int qdir,
    const int qsize,
    const int sid,
    const int sfrom,
    const int chain_score,
    const ChainSeed* chain_seed_array,
    const int chain_seed_count)
{
    int qbeg = chain_seed_array[0].qoff;
    int qend = chain_seed_array[chain_seed_count-1].qoff 
                + 
                chain_seed_array[chain_seed_count-1].length;
    int sbeg = chain_seed_array[0].soff + sfrom;
    int send = chain_seed_array[chain_seed_count-1].soff
                +
                chain_seed_array[chain_seed_count-1].length
                +
                sfrom;
    for (int i = 0; i < hsp_count; ++i) {
        const BlastHSP* hsp = hsp_array + i;
        hbn_assert(hsp->hbn_subject.strand == FWD);
        if (hsp->hbn_query.strand != qdir) continue;
        const int E = 200;
        int r = (qbeg + E >= hsp->hbn_query.offset)
                &&
                (qend <= hsp->hbn_query.end + E)
                &&
                (sbeg + E >= hsp->hbn_subject.offset)
                &&
                (send <= hsp->hbn_subject.end + E);
        if (r) return TRUE;
    }

    return FALSE;
}

static BOOL
hsp_is_contained(const BlastHSP* hsp_array,
    const int hsp_count,
    BlastHSP* newhsp)
{
    const int E = 200;
    for (int i = 0; i < hsp_count; ++i) {
        const BlastHSP* hsp = hsp_array + i;
        if (hsp->hbn_query.strand != newhsp->hbn_query.strand) continue;
        int r = (newhsp->hbn_query.offset + E >= hsp->hbn_query.offset)
                &&
                (newhsp->hbn_query.end <= hsp->hbn_query.end + E)
                &&
                (newhsp->hbn_subject.offset + E >= hsp->hbn_subject.offset)
                &&
                (newhsp->hbn_subject.end <= hsp->hbn_subject.end + E);
        if (r) return TRUE;
    }
    return FALSE;
}

static void s_UpdateEditScript(GapEditScript* esp, int pos, int bf, int af) {
   int op, qd, sd;

   if (bf > 0) {
      op = pos;
      qd = sd = bf;
      do {
          if (--op < 0) return;
          switch(esp->op_type[op]) {
          case eGapAlignSub:
              qd -= esp->num[op];
              sd -= esp->num[op];
              break;
          case eGapAlignIns:
              qd -= esp->num[op];
              break;
          case eGapAlignDel:
              sd -= esp->num[op];
          default:
              break;
          }
      } while (qd > 0 || sd > 0);

      esp->num[op] = -MAX(qd, sd);
      esp->op_type[op++] = eGapAlignSub;
      for (; op < pos-1; op++) esp->num[op] = 0;
      esp->num[pos] += bf;
      qd -= sd;
      esp->op_type[pos-1] = (qd>0) ? eGapAlignDel: eGapAlignIns;
      esp->num[pos-1] = (qd>0) ? qd : -qd;
   }

   if (af > 0) {
      op = pos;
      qd = sd = af;
      do {
          if (++op >= esp->size) return;
          switch(esp->op_type[op]) {
          case eGapAlignSub:
              qd -= esp->num[op];
              sd -= esp->num[op];
              break;
          case eGapAlignIns:
              qd -= esp->num[op];
              break;
          case eGapAlignDel:
              sd -= esp->num[op];
          default:
              break;
          }
      } while (qd > 0 || sd > 0);

      esp->num[op] = -MAX(qd, sd);
      esp->op_type[op--] = eGapAlignSub;
      for (; op > pos+1; op--) esp->num[op] = 0;
      esp->num[pos] += af;
      qd -= sd;
      esp->op_type[pos+1] = (qd>0) ? eGapAlignDel: eGapAlignIns;
      esp->num[pos+1] = (qd>0) ? qd : -qd;
   }
}

static void s_RebuildEditScript(GapEditScript* esp) {
   int i, j;
   for (i=0, j=-1; i<esp->size; i++) {
       if (esp->num[i] == 0) continue;
       if (j>=0 && esp->op_type[i] == esp->op_type[j]) {
           esp->num[j] += esp->num[i];
       } else if (j==-1 || esp->op_type[i] == eGapAlignSub
           || esp->op_type[j] == eGapAlignSub) {
           esp->op_type[++j] = esp->op_type[i];
           esp->num[j] = esp->num[i];
       } else {
           int d = esp->num[j] - esp->num[i];
           if (d > 0) {
              esp->num[j-1] += esp->num[i];
              esp->num[j] = d;
           } else if (d < 0) {
              esp->num[j-1] += esp->num[j];
              esp->num[j] = -d;
              esp->op_type[j] = esp->op_type[i];
           } else {
              esp->num[j-1] += esp->num[j];
              --j;
           }
       }
   }
   esp->size = ++j;
}

static void s_ReduceGaps(GapEditScript* esp, const Uint1 *q, const Uint1 *s,
                                             const Uint1 *qf,const Uint1 *sf){
   int i, j, nm1, nm2, d;
   const Uint1 *q1, *s1;

   for (q1=q, s1=s, i=0; i<esp->size; i++) {
       if (esp->num[i] == 0) continue;
       if (esp->op_type[i] == eGapAlignSub) {
           if(esp->num[i] >= 12) {
               nm1 = 1;
               if (i > 0) {
                   while (q1-nm1>=q && (*(q1-nm1) == *(s1-nm1))) ++nm1;
               }
               q1 += esp->num[i];
               s1 += esp->num[i];
               nm2 = 0;
               if (i < esp->size -1) {
                   while ((q1+1<qf) && (s1+1<sf) && (*(q1++) == *(s1++))) ++nm2;
               }
               if (nm1>1 || nm2>0) s_UpdateEditScript(esp, i, nm1-1, nm2);
               q1--; s1--;
           } else {
               q1 += esp->num[i];
               s1 += esp->num[i];
           }
       } else if (esp->op_type[i] == eGapAlignIns) {
           q1 += esp->num[i];
       } else {
           s1 += esp->num[i];
       }
   }
   s_RebuildEditScript(esp);

   for (i=0; i<esp->size; i++) {
       if (esp->op_type[i] == eGapAlignSub) {
           q += esp->num[i];
           s += esp->num[i];
           continue;
       }
       if (i>1 && esp->op_type[i] != esp->op_type[i-2]
               && esp->num[i-2] > 0) {
           d = esp->num[i] + esp->num[i-1] + esp->num[i-2];
           if (d == 3) {
               /* special case, no need to do further testing */
               (esp->num[i-2]) = 0;
               (esp->num[i-1]) = 2;
               (esp->num[i]) = 0;
               if (esp->op_type[i] == eGapAlignIns) {
                   ++q;
               } else {
                   ++s;
               }
           } else if (d < 12) {
               /* Try reducing this sub... */
               nm1 = 0;
               nm2 = 0;
               d = MIN(esp->num[i], esp->num[i-2]);
               q -= esp->num[i-1];
               s -= esp->num[i-1];
               q1 = q;
               s1 = s;
               if (esp->op_type[i] == eGapAlignIns) {
                   s -= d;
               } else {
                   q -= d;
               }
               for (j=0; j<esp->num[i-1]; ++j, ++q1, ++s1, ++q, ++s) {
                   if (*q1 == *s1) nm1++;
                   if (*q == *s) nm2++;
               }
               for (j=0; j<d; ++j, ++q, ++s) {
                   if (*q == *s) nm2++;
               }
               if (nm2 >= nm1 - d) {
                   (esp->num[i-2]) -= d;
                   (esp->num[i-1]) += d;
                   (esp->num[i]) -= d;
               } else {
                   q = q1;
                   s = s1;
               }
           }
       }
       if (esp->op_type[i] == eGapAlignIns) {
           q += esp->num[i];
       } else {
           s += esp->num[i];
       }
   }
   s_RebuildEditScript(esp);
}

static void 
set_blasthsp(HbnTracebackData* data,
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
                score = (qc == sc) ? reward : penalty;
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
    hsp->query.offset = data->qoff;
    hsp->query.end = data->qend;
    hsp->query.gapped_start = data->qoff;
    hsp->subject.frame = 0;
    hsp->subject.offset = subject_offset + data->soff;
    hsp->subject.end = subject_offset + data->send;
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
    hsp->hbn_query.offset = data->qoff;
    hsp->hbn_query.end = data->qend;
    hsp->hbn_query.seq_size = qsize;

    hsp->hbn_subject.oid = sid;
    hsp->hbn_subject.strand = FWD;
    hsp->hbn_subject.offset = subject_offset + data->soff;
    hsp->hbn_subject.end = subject_offset + data->send;
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
extract_subject_subsequence_without_ambig_res(const text_t* db, const int sid, const size_t from, const size_t to, vec_u8* subject)
{
    kv_clear(*subject);
    const u8* s = db->unpacked_seq + seqdb_seq_offset(db, sid);
    const u8 code_table[16] = {0,1,2,3,0,0,0,0,0,0,0,0,0,0,0,0xf};
    for (size_t i = from; i < to; ++i) {
        u8 c = s[i];
        c = code_table[c];
        kv_push(u8, *subject, c);
    }
}

static void
hbn_extend_subject_subseq_hit_list(HbnSubseqHitExtnData* data,
    const text_t* db,
    HbnSubseqHit* hit_array,
    const int hit_count,
    const int query_id,
    const int fwd_query_context,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_length,
    vec_u8* subject_v,
    const HbnProgramOptions* opts,
    BlastHSPList* hsp_list,
    HbnHSPResults* results)
{
    hsp_list->hbn_best_raw_score = 0;
    hsp_list->hspcnt = 0;
    hsp_list->hsp_array = NULL;
    hsp_list->oid = hit_array[0].sid;
    hsp_list->query_index = query_id;

    sort_subseq_hit_score_gt(hit_count, hit_array);
    MaximalExactMatchWorkData* mem_data = data->mem_data;
    vec_chain_seed* chain_seed_list = &data->chain_seed_list;
    BlastHSP hsp_array[opts->max_hsps_per_subject];
    int hsp_count = 0;
    for (int i = 0; i < hit_count && i < opts->max_hsps_per_subject + 1 && hsp_count < opts->max_hsps_per_subject; ++i) {
        HbnSubseqHit* hit = hit_array + i;
        extract_subject_subsequence_without_ambig_res(db, hit->sid, hit->sfrom, hit->sto, subject_v);
        const u8* subject = kv_data(*subject_v);
        const int subject_length = kv_size(*subject_v);
        int strand = hit->qdir;
        if (!MaximalExactMatchWorkData_FindCandidates(mem_data, subject, subject_length, strand)) {
            continue;
        }
        const int max_k = 5;
        for (size_t k = 0; k < kv_size(mem_data->init_hit_list) && k < max_k; ++k) {
            HbnInitHit* init_hit = kv_data(mem_data->init_hit_list) + k;
            ChainSeed* chain_seed_array = (ChainSeed*)init_hit->chain_seed_array;
            int chain_seed_count = init_hit->chain_seed_count;
            kv_clear(*chain_seed_list);            
            for (int x = 0; x < chain_seed_count; ++x) {
                ChainSeed seed;
                seed.qoff = chain_seed_array[x].soff;
                seed.soff = chain_seed_array[x].qoff;
                seed.length = chain_seed_array[x].length;
                kv_push(ChainSeed, *chain_seed_list, seed);
            }

            if (chain_seed_list_is_contained(hsp_array, 
                    hsp_count, 
                    init_hit->sdir,
                    query_length,
                    hit->sid,
                    hit->sfrom,
                    init_hit->score,
                    kv_data(*chain_seed_list),
                    kv_size(*chain_seed_list))) {
                    continue;
            }
            const u8* query = (init_hit->sdir == FWD) ? fwd_query : rev_query;
            if (!hbn_traceback(data->traceback_data,
                    query,
                    query_length,
                    subject,
                    subject_length,
                    kv_data(*chain_seed_list),
                    kv_size(*chain_seed_list),
                    opts->query_cov_hsp_res,
                    opts->perc_identity,
                    !opts->skip_overhang)) {
                continue;
            }
            BlastHSP* hsp = hsp_array + hsp_count;
            memset(hsp, 0, sizeof(BlastHSP));
            const char* qaln = data->traceback_data->qas;
            const char* saln = data->traceback_data->sas;
            const int aln_len = data->traceback_data->qae - data->traceback_data->qas;
            set_blasthsp(data->traceback_data, query_id, fwd_query_context, init_hit->sdir, query_length, hit->sid, hit->sfrom, seqdb_seq_size(db, hit->sid), qaln, saln, aln_len, opts, hit->score, init_hit->score, hsp);
            s_ReduceGaps(hsp->gap_info,
                query + data->traceback_data->qoff,
                subject + data->traceback_data->soff,
                query + data->traceback_data->qend,
                subject + data->traceback_data->send);
            ++hsp_count;
            if (hsp_count == opts->max_hsps_per_subject) break;
        }
    }
    if (!hsp_count) return;

    ks_introsort_blasthsp_score_gt(hsp_count, hsp_array);
    hsp_list->hbn_best_raw_score = hsp_array[0].score;
    hsp_list->hspcnt = hsp_count;
    hsp_list->hsp_max = hsp_count;
    hsp_list->hsp_array = (BlastHSP**)calloc(hsp_count, sizeof(BlastHSP*));
    for (int i = 0; i < hsp_count; ++i) {
        BlastHSP* hsp = (BlastHSP*)calloc(1, sizeof(BlastHSP));
        memcpy(hsp, hsp_array + i, sizeof(BlastHSP));
        hsp_list->hsp_array[i] = hsp;
    }
}

void
hbn_extend_query_subseq_hit_list(HbnSubseqHit* subseq_hit_array,
    int subseq_hit_count,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_id,
    const int fwd_query_context,
    const int query_length,
    vec_u8* subject_v,
    const text_t* db,
    const HbnProgramOptions* opts,
    HbnSubseqHitExtnData* data,
    BlastHitList* hit_list,
    HbnHSPResults* results)
{
    MaximalExactMatchWorkData_Init(data->mem_data, 
        fwd_query,
        rev_query,
        query_length);
    data->mem_data->ref_name = query_name;

    BlastHSPList hsplist_array[opts->hitlist_size];
    int hsplist_count = 0;
    int i = 0;
    while (i < subseq_hit_count) {
        int j = i + 1;
        while (j < subseq_hit_count && subseq_hit_array[i].sid == subseq_hit_array[j].sid) ++j;
        BlastHSPList* hsp_list = hsplist_array + hsplist_count;
        hbn_extend_subject_subseq_hit_list(data,
            db,
            subseq_hit_array + i,
            j - i,
            query_id,
            fwd_query_context,
            query_name,
            fwd_query,
            rev_query,
            query_length,
            subject_v,
            opts,
            hsp_list,
            results);
        if (hsp_list->hspcnt) ++hsplist_count;
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