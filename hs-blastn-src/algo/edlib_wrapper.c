#include "edlib_wrapper.h"

#include "edlib.h"
#include "hbn_traceback_aux.h"

EdlibAlignData*
EdlibAlignDataNew()
{
    EdlibAlignData* data = (EdlibAlignData*)calloc(1, sizeof(EdlibAlignData));
    ks_init(data->query);
    ks_init(data->target);
    ks_init(data->qaln);
    ks_init(data->taln);
    data->tolerance = -1;
    data->dist = 0;
    data->do_traceback = 1;
    return data;
}

EdlibAlignData*
EdlibAlignDataFree(EdlibAlignData* data)
{
    ks_destroy(data->query);
    ks_destroy(data->target);
    ks_destroy(data->qaln);
    ks_destroy(data->taln);
    free(data);
    return NULL; 
}

int
edlib_nw(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    kstring_t* qaln,
    kstring_t* taln)
{
    ks_clear(*qaln);
    ks_clear(*taln);
    if (query_size == 0 || target_size == 0) return 1;

    ks_set_size(&data->query, query_size);
    for (int i = 0; i < query_size; ++i) {
        int c = query[i];
        hbn_assert(c >= 0 && c < 4);
        ks_A(data->query, i) = DECODE_RESIDUE(c);
    }
    ks_set_size(&data->target, target_size);
    for (int i = 0; i < target_size; ++i) {
        int c = target[i];
        hbn_assert(c >= 0 && c < 4);
        ks_A(data->target, i) = DECODE_RESIDUE(c);
    }
    EdlibAlignTask task = EDLIB_TASK_PATH;
    int tolerance = hbn_max(query_size, target_size) * 0.5;
    EdlibAlignResult align = edlibAlign(ks_s(data->query),
                                query_size,
                                ks_s(data->target),
                                target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task, NULL, 0));

    if (align.numLocations == 0) {
        tolerance = hbn_max(query_size, target_size);
        align = edlibAlign(ks_s(data->query),
                           query_size,
                           ks_s(data->target),
                           target_size,
                           edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task, NULL, 0));        
    }
    hbn_assert(align.numLocations);
    if (align.numLocations == 0) return 0;

    int qBgn = 0;
    int qEnd = query_size;
    int tBgn = align.startLocations[0];
    int tEnd = align.endLocations[0] + 1;
    ks_set_size(qaln, align.alignmentLength + 1);
    ks_set_size(taln, align.alignmentLength + 1);
    edlibAlignmentToStrings(align.alignment,
        align.alignmentLength,
        tBgn,
        tEnd,
        qBgn,
        qEnd,
        ks_s(data->target),
        ks_s(data->query),
        ks_s(*taln),
        ks_s(*qaln));
    edlibFreeAlignResult(align);
    ks_pop_back(*qaln);
    ks_pop_back(*taln);
    hbn_assert(ks_size(*qaln) == ks_size(*taln));
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0,
        query,
        qBgn,
        qEnd,
        ks_s(*qaln),
        0,
        target,
        tBgn,
        tEnd,
        ks_s(*taln),
        ks_size(*qaln),
        TRUE);

    return 1;
}

int
edlib_shw(EdlibAlignData* data,
    const u8* query,
    int query_size,
    const u8* target,
    int target_size,
    int* qend,
    int* tend,
    kstring_t* qaln,
    kstring_t* taln)
{
    ks_clear(*qaln);
    ks_clear(*taln);
    *qend = *tend = 0;
    if (query_size == 0 || target_size == 0) return 0;

    ks_set_size(&data->query, query_size);
    for (int i = 0; i < query_size; ++i) {
        int c = query[i];
        hbn_assert(c >= 0 && c < 4);
        ks_A(data->query, i) = DECODE_RESIDUE(c);
    }
    ks_set_size(&data->target, target_size);
    for (int i = 0; i < target_size; ++i) {
        int c = target[i];
        hbn_assert(c >= 0 && c < 4);
        ks_A(data->target, i) = DECODE_RESIDUE(c);
    }
    EdlibAlignTask task = EDLIB_TASK_PATH;
    int tolerance = hbn_max(query_size, target_size) * 0.35;
    EdlibAlignResult align = edlibAlign(ks_s(data->query),
                                query_size,
                                ks_s(data->target),
                                target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_SHW, task, NULL, 0));
    if (align.numLocations == 0) return 0;

    int qBgn = 0;
    int qEnd = query_size;
    int tBgn = align.startLocations[0];
    int tEnd = align.endLocations[0] + 1;
    hbn_assert(tBgn == 0);
    ks_set_size(qaln, align.alignmentLength + 1);
    ks_set_size(taln, align.alignmentLength + 1);
    edlibAlignmentToStrings(align.alignment,
        align.alignmentLength,
        tBgn,
        tEnd,
        qBgn,
        qEnd,
        ks_s(data->target),
        ks_s(data->query),
        ks_s(*taln),
        ks_s(*qaln));
    edlibFreeAlignResult(align);
    ks_pop_back(*qaln);
    ks_pop_back(*taln);
    hbn_assert(ks_size(*qaln) == ks_size(*taln));
    *qend = qEnd;
    *tend = tEnd;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0,
        query,
        qBgn,
        qEnd,
        ks_s(*qaln),
        0,
        target,
        tBgn,
        tEnd,
        ks_s(*taln),
        ks_size(*qaln),
        TRUE);

    return 1;    
}

static BOOL
get_next_sequence_block(const u8* query,
						int qidx,
						const int qsize,
						const u8* target,
						int tidx,
						const int tsize,
						const int desired_block_size,
						const BOOL right_extend,
						vec_u8* qfrag,
						vec_u8* tfrag)
{
	BOOL last_block = FALSE;
	int qleft = qsize - qidx;
	int tleft = tsize - tidx;
	int qblk;
	int tblk;
	if (qleft < desired_block_size + 100 || tleft < desired_block_size + 100) {
		qblk = tleft * 1.3;
		qblk = hbn_min(qblk, qleft);
		tblk = qleft * 1.3;
		tblk = hbn_min(tblk, tleft);
		last_block = TRUE;
	} else {
		qblk = desired_block_size;
		tblk = desired_block_size;
		last_block = FALSE;
	}
	
	kv_clear(*qfrag);
	kv_clear(*tfrag);
	if (right_extend) {
		const u8* Q = query + qidx;
		for (int i = 0; i < qblk; ++i) kv_push(u8, *qfrag, Q[i]);
		const u8* R = target + tidx;
		for (int i = 0; i < tblk; ++i) kv_push(u8, *tfrag, R[i]);
	} else {
		const u8* Q = query - qidx;
		for (int i = 0; i < qblk; ++i) kv_push(u8, *qfrag, Q[-i]);
		const u8* R = target - tidx;
		for (int i = 0; i < tblk; ++i) kv_push(u8, *tfrag, R[-i]);
	}
	
	return last_block;
}

void
edlib_extend(EdlibAlignData* data,
    const u8* query,
    const int query_size,
    const u8* target,
    const int target_size,
    const int block_size,
    const BOOL right_extend,
    vec_u8* qfrag,
    vec_u8* tfrag,
    int* qend,
    int* tend,
    kstring_t* qaln,
    kstring_t* taln)
{
    ks_clear(*qaln);
    ks_clear(*taln);
    int qidx = 0, tidx = 0;
    while (1) {
        int qfae, tfae, qfrag_size, tfrag_size;
        BOOL last_block = get_next_sequence_block(query,
                            qidx,
                            query_size,
                            target,
                            tidx,
                            target_size,
                            block_size,
                            right_extend,
                            qfrag,
                            tfrag);
        qfrag_size = kv_size(*qfrag);
        tfrag_size = kv_size(*tfrag);
        if (qfrag_size == 0 || tfrag_size == 0) break;

        edlib_shw(data,
            kv_data(*qfrag),
            qfrag_size,
            kv_data(*tfrag),
            tfrag_size,
            &qfae,
            &tfae,
            &data->qaln,
            &data->taln);

        BOOL done = last_block;
        if (qfrag_size - qfae > 30 || tfrag_size - tfae > 30) done = TRUE;
        int acnt = 0, qcnt = 0, tcnt = 0;
        const int M = 8;
        int align_size = ks_size(data->qaln);
        int k = align_size - 1, m = 0;
        while (k >= 0) {
            const char qc = ks_A(data->qaln, k);;
            const char tc = ks_A(data->taln, k);
            if (qc != GAP_CHAR) ++qcnt;
            if (tc != GAP_CHAR) ++tcnt;
            m = (qc == tc) ? (m+1) : 0;
            ++acnt;
            if (m == M) break;
            --k;
        }

        if (m != M || k < 1) {
            align_size = 0;
            ks_clear(data->qaln);
            ks_clear(data->taln);
            for (int i = 0; i < qfrag_size && i < tfrag_size; ++i) {
                int qc = kv_A(*qfrag, i);
                int tc = kv_A(*tfrag, i);
                if (qc != tc) break;
                qc = DECODE_RESIDUE(qc);
                tc = DECODE_RESIDUE(tc);
                kputc(qc, &data->qaln);
                kputc(tc, &data->taln);
                ++align_size;
            }
            done = TRUE;
        } else {
            align_size -= acnt;
            qidx += (qfae - qcnt);
            tidx += (tfae - tcnt);
            if (done) align_size += M;
        }

        kputsn(ks_s(data->qaln), align_size, qaln);
        kputsn(ks_s(data->taln), align_size, taln);
        if (done) break;
    }

	int qe = 0, te = 0;
	for (size_t i = 0; i != ks_size(*qaln); ++i) {
		if (ks_A(*qaln, i) != GAP_CHAR) ++qe;
		if (ks_A(*taln, i) != GAP_CHAR) ++te;
	}
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0,
        query,
        0,
        qe,
        ks_s(*qaln),
        0,
        target,
        0,
        te,
        ks_s(*taln),
        ks_size(*qaln),
        right_extend);
    *qend = qe;
    *tend = te;
}