#ifndef GAPALIGN_H
#define	GAPALIGN_H

#include "def.h"
#include "greedy_align.h"
#include "memallocator.h"

enum EGapAlignOpType {
    eGapAlignDel = 0, /**< Deletion: a gap in query */
    eGapAlignSub = 1, /**< Substitution */
    eGapAlignIns = 2, /**< Insertion: a gap in subject */
    eGapAlignInvalid = 3 /**< Invalid operation */
};

struct GapEditScript {
    EGapAlignOpType* op_type; /**< Array of type of operation */
    Int4* num; /**< Array of number of operations */
    Int4 size; /**< Size of above arrays. */
    Int4 num_alloc;
    void Destroy();
};

GapEditScript*
GapEditScriptNew(SmallObjAllocator& soa, Int4 size);

GapEditScript*
GapEditScriptDelete(SmallObjAllocator& soa, GapEditScript* esp);

void GapEditScriptPrint(GapEditScript* esp);

Int2 GapEditScriptPartialCopy(GapEditScript* new_esp, 
                              int offset, 
                              const GapEditScript* old_esp,
                              int start,
                              int stop);

struct GapPrelimEditScript {
    EGapAlignOpType op_type; /**< Type of operation */
    Int4 num; /**< Number of operations */
};

struct GapPrelimEditBlock {
    GapPrelimEditScript *edit_ops; /**< array of edit operations */
    Int4 num_ops_allocated; /**< size of allocated array */
    Int4 num_ops; /**< number of edit ops presently in use */
    EGapAlignOpType last_op; /**< most recent operation added */

    Int2 Realloc(Int4 total_ops);
    Int2 AddNew(EGapAlignOpType op_type, Int4 num);
    void Add(EGapAlignOpType op_type, Int4 num);

    GapPrelimEditBlock();
    void Destroy();
    ~GapPrelimEditBlock();
    void Reset();
    
    void Print();
};

/** Structure supporting the gapped alignment */
// See BlastGapAlignStruct
struct GapAlignStruct {
    GapEditScript* edit_srcipt;
    GapPrelimEditBlock* fwd_prelim_tback;
    GapPrelimEditBlock* rev_prelim_tback;
    SGreedyAlignMem* greedy_align_mem;

    Int4 query_start;
    Int4 query_stop;
    Int4 subject_start;
    Int4 subject_stop;
    Int4 greedy_query_seed_start;
    Int4 greedy_subject_seed_start;
    Int4 score;

    BlastScoringParameters* bsp;
    BlastExtensionParameters* bep;

    GapAlignStruct(BlastScoringParameters* score_params,
            BlastExtensionParameters* ext_params);
    void GapAlignStructFree();
    ~GapAlignStruct();
};

struct HSP {
    Int4 q_off, q_end;
    Int8 s_off, s_end;

    Int4 subject_id;
    Int4 context;
    
    Int4 query_gapped_start;
    Int8 subject_gapped_start;

    Int4 score;
    Int4 num_ident;
    double bit_score;
    double evalue;

    GapEditScript* esp;
    
    Int2 s_Blast_HSPGetNumIdentitesAndPositives(const Uint1* query, const Uint1* subject,
                Int4* num_ident_ptr, Int4* align_length_ptr, Int4* num_pos_ptr);
};

HSP* HSPNew(SmallObjAllocator& soa);
HSP* HSPDelete(SmallObjAllocator& soa, HSP* hsp);
void HSPPrint(HSP* hsp);

/** Cleans out the NULLed out HSP's from the HSP array that
 * is part of the BlastHSPList.
 * @param hsp_list Contains array of pointers to HSP structures [in]
 * @return status of function call.
*/
Int2 Blast_HSPListPurgeNullHSPs(HSP** hsp_list, Int4& num_hsps);

/** Check for an overlap of two different alignments and remove redundant HSPs.
 * A sufficient overlap is when two alignments have the same start or end values
 * If an overlap is found the HSP with the lowest score is removed, if both scores
 * are the same then the first is removed.
 * @param program Type of BLAST program. For some programs (PHI BLAST), the
 *                purge should not be performed. [in]
 * @param hsp_list Contains array of pointers to HSPs to purge [in]
 * @param purge Should the hsp be purged? [in]
 * @return The number of valid alignments remaining. 
*/
Int4
Blast_HSPListPurgeHSPsWithCommonEndpoints(HSP** hsp_list,
                                          Int4& num_hsps,
                                          Boolean purge,
                                          SmallObjAllocator& soa);

int ScoreCompareHSPs(const void* h1, const void* h2);
int EvalueCompareHSPs(const void* v1, const void* v2);

class GreedyAligner {
public:
    GreedyAligner(BlastScoringParameters* score_params,
            BlastExtensionParameters* ext_params,
            SmallObjAllocator& alloc);
    ~GreedyAligner();
    Int2 GreedyGappedAlignment(const Uint1* query, const Uint1* subject, Int4 query_length,
            Int4 subject_length, Int4 q_off, Int4 s_off, Boolean do_traceback, int mark);
    void PackHSP(HSP& hsp, Int4 q_off, Int8 s_off);

private:

    void PrelimEditBlockToGapEditScript(GapPrelimEditBlock* rev_prelim_tback,
            GapPrelimEditBlock* fwd_prelim_tback);

    Int2 GreedyGapAlignStructFill(Int4 q_start, Int4 s_start,
            Int4 q_end, Int4 s_end,
            Int4 q_seed_start,
            Int4 s_seed_start, Int4 score);
    
public:
    GapAlignStruct gap_align;
    GapEditScript* esp;
    SmallObjAllocator& soa;
    static const Int4 kInvalidOffset = -2;
};

#endif	/* GAPALIGN_H */

