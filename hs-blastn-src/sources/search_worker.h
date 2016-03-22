#ifndef SEARCH_WORKER_H
#define SEARCH_WORKER_H

#include "index.h"
#include "symdust.hpp"
#include "arguments.h"
#include "query_info.h"
#include "parameters.h"
#include "gapalign.h"
#include "interval_tree.h"
#include "traceback.h"
#include "memallocator.h"
#include "result_format.h"
#include "seq_masker.hpp"

// Store the information for one ungapped alignment
struct Hit
{
    Int4 soff;
    Int4 qoff;
    Int4 len;
    Int4 score;
    Int4 context;
};

struct SearchWorker
{
	// Memory pool
    SmallObjAllocator soa;
    QueryInfo* global_queries;
	QueryInfo local_queries;
    FMIndex* fmindex;
    DbInfo* dbinfo;
    CSymDustMasker* dust_masker;
    CSeqMasker* window_masker;
	OutputFormat* results;
    Options* options;
    BlastScoreBlk* sbp;
    BlastEffectiveLengthsOptions* eff_options;
	int thread_id;

	// various parameters (used in NCBI-BLAST)
    BlastInitialWordParameters* word_params;
    BlastExtensionParameters* ext_params;
    BlastHitSavingParameters* hit_params;
    BlastEffectiveLengthsParameters* eff_params;
    BlastScoringParameters* score_params;

    SearchWorker(Options* opts, 
                 QueryInfo* query_info_,
                 FMIndex* index,
                 DbInfo* di,
                 CSymDustMasker* csdm,
                 CSeqMasker* win_m,
				 int tid);
    ~SearchWorker();

	/// Step 1
	// setup the information for each query.
	// e.g., the x-drop, effective search space,...
    Int2 BLAST_GapAlignSetup();
    
	/// Step 2
	// seeding
    void Seeding();
    
    cy_utility::SimpleArray<MEM> seeds;
    cy_utility::SimpleArray<SearchInterval> sis;    
    
	/// Step 3
	// preliminary search 
    cy_utility::SimpleArray<Hit> ungapped_results;
    cy_utility::SimpleArray<HSP*> gapped_results;
    GreedyAligner* gapped_aligner;
	// (two-dimensional) interval tree
    IntervalTree itree;
    void GetGappedScore();
    
	/// Step 4
	// traceback search
    void BLAST_ComputeTraceback();
    
	/// Step 5
	// Output results
    void OutputResult(OutputFormat* result_writter);
    void CleanUp();
    
	// Steps 1-5 are performed in order in this function
    void Go();
	
	void FlushResults(FILE* file);
};

void PrelimSearchStage(cy_utility::SimpleArray<MEM>& seeds,
					   QueryInfo* query_info,
					   DbInfo* dbinfo,
					   Int4 seed_size,
					   cy_utility::SimpleArray<Hit>& ungapped_alignments,
					   cy_utility::SimpleArray<HSP*>& gapped_alignments,
					   BlastInitialWordParameters* word_params,
					   BlastScoreBlk* sbp,
					   BlastHitSavingParameters* hit_params,
					   IntervalTree* itree,
					   Int4 min_diag_separation,
					   GreedyAligner* gapped_aligner,
					   SmallObjAllocator& soa);


#endif // SEARCH_WORKER_H
