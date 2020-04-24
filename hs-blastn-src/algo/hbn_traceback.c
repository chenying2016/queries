#include "hbn_traceback.h"

#include <math.h>

#include "hbn_traceback_aux.h"
#include "mem_finder.h"

static const int kMatLen = 8;
//static const int kMaxEdlibOverHang = 3000;
static const int kMaxOverHang = 1000;

HbnTracebackData*
HbnTracebackDataNew()
{
    HbnTracebackData* data = (HbnTracebackData*)calloc(1, sizeof(HbnTracebackData));
    ks_init(data->qabuf);
    ks_init(data->sabuf);
    ks_init(data->ext_qabuf);
    ks_init(data->ext_sabuf);
    kv_init(data->qfrag);
    kv_init(data->sfrag);
    kv_init(data->trace_seeds);
    data->edlib = EdlibAlignDataNew();
    data->dalign = DalignDataNew(0.35);
    return data;
}

HbnTracebackData*
HbnTracebackDataFree(HbnTracebackData* data)
{
    ks_destroy(data->qabuf);
    ks_destroy(data->sabuf);
    ks_destroy(data->ext_qabuf);
    ks_destroy(data->ext_sabuf);
    kv_destroy(data->qfrag);
    kv_destroy(data->sfrag);
    kv_destroy(data->trace_seeds);
    EdlibAlignDataFree(data->edlib);
    DalignDataFree(data->dalign);
    free(data);
    return NULL;
}

static void
compute_trace_points(const ChainSeed* seed_array,
    const int seed_count,
    vec_chain_seed* trace_seeds)
{
    kv_clear(*trace_seeds);
    if (seed_count == 0) return;
    int i = 0;
    while (i < seed_count) {
        int f = i;
        while (f < seed_count - 1) {
            int g = f + 1;
            ChainSeed sf = seed_array[f];
            ChainSeed sg = seed_array[g];
            hbn_assert(sf.qoff < sg.qoff, 
                "f = %d, g = %d, sf = (%d, %d, %d), sg = (%d, %d, %d)", 
                f, g, sf.qoff, sf.soff, sf.length, sg.qoff, sg.soff, sg.length);
            hbn_assert(sf.soff < sg.soff,
                "seed_count = %d, f = %d, g = %d, sf = (%d, %d, %d), sg = (%d, %d, %d)", 
                seed_count, f, g, sf.qoff, sf.soff, sf.length, sg.qoff, sg.soff, sg.length);
            int dq = sg.qoff - sf.qoff;
            int ds = sg.soff - sf.soff;
            int r1 = abs(dq - ds) <= 100;
            int r2 = fabs(1.0 - 1.0 * dq / ds) <= 0.21;
            if ((!r1) || (!r2)) break;
            ++f;            
        }

        if (f > i) {
            for (int k = i; k <= f; ++k) {
                kv_push(ChainSeed, *trace_seeds, seed_array[k]);
            }
        } else if (i == 0) {
            kv_push(ChainSeed, *trace_seeds, seed_array[0]);
        }

        i = f + 1;
    }

    if (kv_back(*trace_seeds).qoff < seed_array[seed_count-1].qoff) {
        kv_push(ChainSeed, *trace_seeds, seed_array[seed_count-1]);
    }
}

static void
HbnTracebackDataInit(HbnTracebackData* data,
    int qoff,
    int qsize,
    int soff,
    int ssize)
{
    int wrk_l = 2 * hbn_max(qsize, ssize);
    ks_reserve(&data->qabuf, wrk_l);
    ks_reserve(&data->sabuf, wrk_l);
    int wrk_ll = 2 * hbn_max(qoff, soff);
    data->qas = data->qae = ks_s(data->qabuf) + wrk_ll;
    data->sas = data->sae = ks_s(data->sabuf) + wrk_ll;
    *data->qae = '\0';
    *data->sae = '\0';
    data->qoff = data->qend = 0;
    data->soff = data->send = 0;
    data->dist = 0;
    data->ident_perc = 0;
}

static void
apped_match_subseq(const u8* query,
    const int qfrom,
    const int qto,
    const u8* subject,
    const int sfrom,
    const int sto,
    char** qae,
    char** sae)
{
    //HBN_LOG("qf = %d, qt = %d, sf = %d, st = %d", qfrom, qto, sfrom, sto);
    hbn_assert(qto - qfrom == sto - sfrom);
    int n = qto - qfrom;
    for (int i = 0; i < n; ++i) {
        int c1 = query[qfrom + i];
        c1 = DECODE_RESIDUE(c1);
        **qae = c1;
        ++(*qae);
        int c2 = subject[sfrom + i];
        c2 = DECODE_RESIDUE(c2);
        **sae = c2;
        ++(*sae);
        hbn_assert(c1 == c2, "qfrom = %d, qto = %d, sfrom = %d, sto = %d, qi = %d, si = %d, c1 = %c, c2 = %c", 
            qfrom, qto, sfrom, sto,
            qfrom + i, sfrom + i, c1, c2);
    }
    **qae = '\0';
    **sae = '\0';
}

static int
run_nw(const u8* query, 
    const int query_length,
    const u8* subject, 
    const int subject_length,
    HbnTracebackData* data,
    kstring_t* qaln,
    kstring_t* saln,
    char** qae,
    char** sae)
{
    edlib_nw(data->edlib, query, query_length, subject, subject_length, qaln, saln);
    hbn_assert(ks_size(*qaln) == ks_size(*saln));
    memcpy(*qae, ks_s(*qaln), ks_size(*qaln));
    *qae += ks_size(*qaln);
    **qae = '\0';
    memcpy(*sae, ks_s(*saln), ks_size(*saln));
    *sae += ks_size(*saln);
    **sae = '\0';
    return 1;
}

int
overhang_extend(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    kstring_t* qaln,
    kstring_t* saln)
{
    int qoff, qend, soff, send;
    double ident_perc;
    int r = dalign_align(dalign,
                query,
                0,
                query_length,
                subject,
                0,
                subject_length,
                1,
                0.65,
                &qoff,
                &qend,
                &soff,
                &send,
                &ident_perc,
                NULL,
                NULL);
    if (!r) return r;
    hbn_assert(qoff == 0 && soff == 0);

    r = edlib_nw(edlib,
            query + qoff,
            qend - qoff,
            subject + soff,
            send - soff,
            qaln,
            saln);
    return r;
}

static int
left_extend(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    const u8* subject,
    vec_u8* qsbuf,
    vec_u8* ssbuf,
    kstring_t* qabuf,
    kstring_t* sabuf,
    int* qbeg,
    int* sbeg,
    char** qas,
    char** sas)
{
    int qls = *qbeg;
    int sls = *sbeg;
    int ls = hbn_min(qls, sls);
    if (ls > kMaxOverHang) return 0;
    if (ls == 0) return 0;

    qls += kMatLen;
    sls += kMatLen;

    kv_clear(*qsbuf);
    kv_clear(*ssbuf);
    for (ls = qls; ls; --ls) kv_push(u8, *qsbuf, query[ls-1]);
    for (ls = sls; ls; --ls) kv_push(u8, *ssbuf, subject[ls-1]);

    int r = overhang_extend(dalign, edlib, kv_data(*qsbuf), qls, kv_data(*ssbuf), sls, qabuf, sabuf);
    if (!r) return r;
    *qas += kMatLen;
    *sas += kMatLen;
    *qbeg += kMatLen;
    *sbeg += kMatLen;
    int qi = 0, si = 0;
    hbn_assert(ks_size(*qabuf) == ks_size(*sabuf));
    for (size_t i = 0; i < ks_size(*qabuf); ++i) {
        char qc = ks_A(*qabuf, i);
        --(*qas);
        **qas = qc;
        if (qc != GAP_CHAR) ++qi;
        char sc = ks_A(*sabuf, i);
        --(*sas);
        **sas = sc;
        if (sc != GAP_CHAR) ++si;
    }

    *qbeg -= qi;
    *sbeg -= si;
    return 1;
}

static int
right_extend(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    int* qend,
    const int query_length,
    const u8* subject,
    int* send,
    const int subject_length,
    kstring_t* qabuf,
    kstring_t* sabuf,
    char** qae,
    char** sae)
{
    int qrs = query_length - (*qend);
    int srs = subject_length - (*send);
    int rs = hbn_min(qrs, srs);
    if (rs > kMaxOverHang || rs == 0) return 0;

    qrs += kMatLen;
    srs += kMatLen;
    const u8* q = query + query_length - qrs;
    const u8* s = subject + subject_length - srs;

    int r = overhang_extend(dalign, edlib, q, qrs, s, srs, qabuf, sabuf);
    if (!r) return r;

    *qend -= kMatLen;
    *send -= kMatLen;
    *qae -= kMatLen;
    *sae -= kMatLen;
    int qi = 0;
    int si = 0;
    hbn_assert(ks_size(*qabuf) == ks_size(*sabuf));
    for (size_t i = 0; i < ks_size(*qabuf); ++i) {
        char qc = ks_A(*qabuf, i);
        **qae = qc; ++(*qae);
        if (qc != GAP_CHAR) ++qi;
        char sc = ks_A(*sabuf, i);
        **sae = sc; ++(*sae);
        if (sc != GAP_CHAR) ++si;
    }
    **qae = '\0';
    **sae = '\0';
    *qend += qi;
    *send += si;
    return 1;
}

int
hbn_traceback(HbnTracebackData* data,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    const ChainSeed* seed_array,
    const int seed_count,
    const int min_align_size,
    const double min_ident_perc,
    const int process_over_hang)
{
    validate_mem(HBN_LOG_ARGS_DEFAULT, query, subject, seed_array, seed_count);
    const int kTracebackExtendBlock = 1024;
    compute_trace_points(seed_array, seed_count, &data->trace_seeds);
    ChainSeed* tsa = kv_data(data->trace_seeds);
    int tsc = kv_size(data->trace_seeds);
    HbnTracebackDataInit(data, tsa[0].qoff, query_length, tsa[0].soff, subject_length);
    char** qae = &data->qae;
    char** sae = &data->sae;
    const int E = 10;
    const int E2 = E * 2;
    int qfrom = 0, qto = 0;
    int sfrom = 0, sto = 0;
    int r = 1;
    if (tsa[0].length > E2) {
        qfrom = tsa[0].qoff;
        qto = tsa[0].qoff + tsa[0].length;
        sfrom = tsa[0].soff;
        sto = tsa[0].soff + tsa[0].length;
        apped_match_subseq(query, qfrom, qto - E, subject, sfrom, sto - E, qae, sae);
        qfrom = qto - E;
        sfrom = sto - E;
    } else {
        qfrom = tsa[0].qoff;
        sfrom = tsa[0].soff;
    }

    for (int i = 0; i < tsc - 1; ++i) {
        ChainSeed sj = tsa[i+1];
        //HBN_LOG("%d: (%d, %d, %d)", i, sj.qoff, sj.soff, sj.length);
        if (sj.length > E2) {
            qto = sj.qoff + E;
            sto = sj.soff + E;
            run_nw(query + qfrom, qto - qfrom, subject + sfrom, sto - sfrom, data, &data->ext_qabuf, &data->ext_sabuf, qae, sae);
            hbn_assert(r);
            qfrom = qto;
            qto = sj.qoff + sj.length - E;
            sfrom = sto;
            sto = sj.soff + sj.length - E;
            apped_match_subseq(query, qfrom, qto, subject, sfrom, sto, qae, sae);
            qfrom = qto;
            sfrom = sto;
        } else {
            qto = sj.qoff + sj.length / 2;
            sto = sj.soff + sj.length / 2;
            run_nw(query + qfrom, qto - qfrom, subject + sfrom, sto - sfrom, data, &data->ext_qabuf, &data->ext_sabuf, qae, sae);
            hbn_assert(r);
            qfrom = qto;
            sfrom = sto;
        }
    }

    ChainSeed se = tsa[tsc-1];
    qto = se.qoff + se.length;
    sto = se.soff + se.length;
    apped_match_subseq(query, qfrom, qto, subject, sfrom, sto, qae, sae);
    hbn_assert(strlen(data->qas) == strlen(data->sas));

    data->qoff = tsa[0].qoff;
    data->qend = qto;
    data->qsize = query_length;
    data->soff = tsa[0].soff;
    data->send = sto;
    data->ssize = subject_length;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);
    
    int qcnt = 0, scnt = 0;
    char* qa = data->qas;
    char* sa = data->sas;
    int m = 0;
    while (qa < data->qae) {
        int qc = *qa; ++qa;
        int sc = *sa; ++sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    if (m < kMatLen) return 0;
    data->qas = qa - kMatLen;
    data->sas = sa - kMatLen;
    data->qoff += qcnt - kMatLen;
    data->soff += scnt - kMatLen;   

    qcnt = 0;
    scnt = 0;
    qa = data->qae;
    sa = data->sae;
    m = 0;
    while (qa > data->qas) {
        --qa;
        --sa;
        int qc = *qa;
        int sc = *sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    hbn_assert(m == kMatLen, "m = %d, qcnt = %d, scnt = %d", m, qcnt, scnt);
    if (qa == data->qas) return 0;

    data->qae = qa + kMatLen;
    data->sae = sa + kMatLen;
    data->qend -= qcnt - kMatLen;
    data->send -= scnt - kMatLen;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);

    if (1) { //(data->qoff <= kMaxEdlibOverHang || data->soff <= kMaxEdlibOverHang) {
    edlib_extend(data->edlib,
        query + data->qoff - 1,
        data->qoff,
        subject + data->soff - 1,
        data->soff,
        kTracebackExtendBlock,
        FALSE,
        &data->qfrag,
        &data->sfrag,
        &qcnt,
        &scnt,
        &data->ext_qabuf,
        &data->ext_sabuf);

    hbn_assert(qcnt <= data->qoff);
    hbn_assert(scnt <= data->soff);
    hbn_assert(ks_size(data->ext_qabuf) == ks_size(data->ext_sabuf));
    for (size_t i = 0; i < ks_size(data->ext_qabuf); ++i) {
        int res = ks_A(data->ext_qabuf, i);
        --data->qas;
        *data->qas = res;
        res = ks_A(data->ext_sabuf, i);
        --data->sas;
        *data->sas = res;
    }
    data->qoff -= qcnt;
    data->soff -= scnt;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);
    }

    if (1) { //if (query_length - data->qend <= kMaxEdlibOverHang || subject_length - data->send <= kMaxEdlibOverHang) {
    edlib_extend(data->edlib,
        query + data->qend,
        query_length - data->qend,
        subject + data->send,
        subject_length - data->send,
        kTracebackExtendBlock,
        TRUE,
        &data->qfrag,
        &data->sfrag,
        &qcnt,
        &scnt,
        &data->ext_qabuf,
        &data->ext_sabuf);
    hbn_assert(data->qend + qcnt <= query_length);
    hbn_assert(data->send + scnt <= subject_length);
    hbn_assert(ks_size(data->ext_qabuf) == ks_size(data->ext_sabuf));
    memcpy(data->qae, ks_s(data->ext_qabuf), ks_size(data->ext_qabuf));
    data->qae += ks_size(data->ext_qabuf);
    *data->qae = '\0';
    memcpy(data->sae, ks_s(data->ext_sabuf), ks_size(data->ext_sabuf));
    data->sae += ks_size(data->ext_sabuf);
    *data->sae = '\0';
    data->qend += qcnt;
    data->send += scnt;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);
    }

    if (data->qas == data->qae) return 0;

    if (process_over_hang) {
        //HBN_LOG("before:");
        //HbnTracebackDataDump(fprintf, stderr, data);
        //left_extend(data->ksw, query, subject, &data->ksw->qfrag, &data->ksw->tfrag,
        //    &data->qoff, &data->soff, &data->qas, &data->sas);
        //right_extend(data->ksw, query, &data->qend, query_length,
        //    subject, &data->send, subject_length, &data->qae, &data->sae);

        left_extend(data->dalign, data->edlib, query, subject, &data->qfrag, &data->sfrag,
            &data->ext_qabuf, &data->ext_sabuf, &data->qoff, &data->soff, &data->qas, &data->sas);
        right_extend(data->dalign, data->edlib, query, &data->qend, query_length,
            subject, &data->send, subject_length, &data->ext_qabuf, &data->ext_sabuf,
            &data->qae, &data->sae);
        //HBN_LOG("after:");
        //HbnTracebackDataDump(fprintf, stderr, data);
    }

    data->ident_perc = calc_ident_perc(data->qas, data->sas, strlen(data->qas), &data->dist, &data->score);

    truncate_align_bad_ends(data->qas,
        data->sas,
        data->qae - data->qas,
        1,
        &data->qoff,
        &data->qend,
        &data->soff,
        &data->send,
        &data->qas,
        &data->qae,
        &data->sas,
        &data->sae);

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);
    
    r = (data->qae - data->qas >= min_align_size) && (data->ident_perc >= min_ident_perc);
    return r;
}

BOOL
truncate_align_bad_ends(const char* qaln,
    const char* saln,
    const int aln_size,
    int mat_len,
    int* qoff,
    int* qend,
    int* soff,
    int* send,
    const char** qas_,
    const char** qae_,
    const char** sas_,
    const char** sae_)
{
    const char* const qas = qaln;
    const char* const qae = qaln + aln_size;
    const char* const sas = saln;
    const char* const sae = saln + aln_size;
    int qcnt = 0, scnt = 0;
    const char* qa = qas;
    const char* sa = sas;
    int m = 0;

    while (qa < qae) {
        int qc = *qa; ++qa;
        int sc = *sa; ++sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == mat_len) break;
    }
    if (m < mat_len || qa == qae) return 0;
    *qas_ = qa - mat_len;
    *sas_ = sa - mat_len;
    *qoff += qcnt - mat_len;
    *soff += scnt - mat_len;  

    qcnt = 0;
    scnt = 0;
    qa = qae;
    sa = sae;
    m = 0;
    while (qa > qas) {
        --qa;
        --sa;
        int qc = *qa;
        int sc = *sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == mat_len) break;
    }
    hbn_assert(m == mat_len, "m = %d, qcnt = %d, scnt = %d", m, qcnt, scnt);
    if (qa == qas) return 0;

    *qae_ = qa + mat_len;
    *sae_ = sa + mat_len;
    *qend -= qcnt - mat_len;
    *send -= scnt - mat_len; 
    return TRUE;
}