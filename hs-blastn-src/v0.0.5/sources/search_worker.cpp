#include "search_worker.h"
#include "result_format.h"
#include "utility.h"

SearchWorker::SearchWorker(Options* opts, 
                           QueryInfo* qi,
                           FMIndex* index,
                           DbInfo* di,
                           CSymDustMasker* csdm,
                           CSeqMasker* win_m,
						   int tid)
{
    ASSERT(opts != NULL);
    global_queries = qi;
    options = opts;
    fmindex = index;
    dbinfo = di;
    results = new OutputFormat(options->output_options, options->hit_options, dbinfo);
    dust_masker = csdm;
    window_masker = win_m;
	thread_id = tid;
    
    word_params = NULL;
    ext_params = NULL;
    hit_params = NULL;
    eff_params = NULL;
    score_params = NULL;
    
    sbp = NULL;
    gapped_aligner = NULL;
    
    eff_options = BlastEffectiveLengthsOptionsNew();
}

SearchWorker::~SearchWorker()
{
    if (sbp) delete sbp;

    if (word_params) delete word_params;
    if (ext_params) delete ext_params;
    if (hit_params) delete hit_params;
    if (eff_params) delete eff_params;
    if (score_params) delete score_params;   
    
    if (gapped_aligner) delete gapped_aligner;
    
    soa.Release();
    
    if (eff_options) free(eff_options);
	if (results) delete results;
}

Int2 SearchWorker::BLAST_GapAlignSetup()
{
    Int2 status = 0;
    Int8 max_subject_length = 0;
    Int8 min_subject_length = 0;
    Int8 total_length = -1;
    Int4 num_seqs = -1;
    
    sbp = new BlastScoreBlk(BLASTNA_SEQ_CODE, &local_queries);
    sbp->ScoreBlkInit(options->scoring_options, 1.0);
    
    ASSERT(dbinfo != NULL);
    total_length = dbinfo->GetDbLength();
    ASSERT(total_length > 0);
    num_seqs = dbinfo->GetNumSeqs();
    ASSERT(num_seqs > 0);

    eff_params = new BlastEffectiveLengthsParameters(eff_options, total_length, num_seqs);
    status = sbp->BLAST_CalcEffLength(options->scoring_options, eff_params, &local_queries);
    ASSERT(status == 0);

    score_params = new BlastScoringParameters(options->scoring_options, sbp);
    ext_params = new BlastExtensionParameters(options->ext_options, sbp, &local_queries);
    
    min_subject_length = total_length / num_seqs;
    hit_params = new BlastHitSavingParameters(options->hit_options, sbp, &local_queries, min_subject_length);
    
    word_params = new BlastInitialWordParameters(options->word_options, hit_params, sbp,
                                                 &local_queries, total_length);
    
    eff_options->num_searchspaces = sbp->number_of_contexts;
    eff_options->searchsp_eff = sbp->eff_searchsp;
    
    gapped_aligner = new GreedyAligner(score_params, ext_params, soa);
    
    return status;
}

void SearchWorker::Seeding()
{
    Uint1* q;
    Int4 qlen;
    Uint4 n_ranges = local_queries.scan_ranges.size();;
    Uint4 i;
    Int4 context = 0;
    
    sis.clear();

    for (i = 0; i < n_ranges; ++i)
    {
        Point2d_Int4& range = local_queries.scan_ranges[i];
        q = (Uint1*)&local_queries.blastna_query[range[0]];
        qlen = range[1] - range[0] + 1;
        Int8 range_start = static_cast<Int8>(range[0]);
        while (range_start > local_queries.contexts[context].offset + local_queries.contexts[context].length)++context;
        
        Int4 context_start = local_queries.GetContextOffset(context);
        
        Int4 start_scan = 0, end_scan;
		while (start_scan < qlen)
		{
			while (start_scan < qlen && q[start_scan] > 3) 
				++start_scan;
			end_scan = start_scan + 1;
			while (end_scan < qlen && q[end_scan] < 4)
				++end_scan;
			fmindex->Seeding(q + start_scan, end_scan - start_scan,
                         context, options->seed_options->seed_size,
                         range[0] + start_scan - context_start, sis);
			start_scan = end_scan;
		}
    } 
    
    seeds.clear();
	int n;
    n = fmindex->LocateSeeds(sis, seeds, dbinfo->GetDbLength(), 
                         options->seed_options->seed_size,
                         &local_queries);

	cy_utility::Log::Trace(__func__, "Number of seeds: %d", seeds.size());
}

void SearchWorker::BLAST_ComputeTraceback()
{
    Int4 i, j, k;
    Int4 context;
    if (gapped_results.size() == 0) return;
    
    i = 0;
    Int4 num_hsps = gapped_results.size();
    HSP** context_hsps;
    Int4 num_context_hsps;
    
    Int4 next_final_gapped_result_index = 0;
    HSP** final_gapped_results = (HSP**)gapped_results.get_data();
    
    while (i < num_hsps)
    {
        ASSERT(gapped_results[i] != NULL);
        context = gapped_results[i]->context;  
        num_context_hsps = gapped_results[i]->num_ident;
        context_hsps = &gapped_results[i];
        i += num_context_hsps;

        Blast_TraceBackFromOneContextHSPList(context_hsps, 
                                             num_context_hsps,
                                             soa,
                                             &local_queries,
                                             sbp,
                                             dbinfo,
                                             itree,
                                             gapped_aligner,
                                             score_params,
                                             ext_params,
                                             hit_params);          

        for (k = 0; k < num_context_hsps; ++k)
        {
            final_gapped_results[next_final_gapped_result_index + k]
                    = context_hsps[k];
        }
        if (num_context_hsps > 0)
        {
            final_gapped_results[next_final_gapped_result_index]->query_gapped_start
                    = num_context_hsps;
        }
        next_final_gapped_result_index += num_context_hsps;
    }

    gapped_results.set_array_size(next_final_gapped_result_index);

	cy_utility::Log::Trace(__func__, "Final gapped alignments: %d", next_final_gapped_result_index);
}

void SearchWorker::OutputResult(OutputFormat* result_writter)
{
    Int4 qid;
    Int4 i = 0, j, k;
    Int4 num_hsps = gapped_results.size();
    HSP* hsp;  

    while (i < num_hsps)
    {
        qid = (gapped_results[i]->context) / 2;
		local_queries.query_offsets[qid].result_offset = i;
        j = gapped_results[i]->query_gapped_start;
        HSP** hsps = (HSP**)gapped_results.get_data() + i;
        
        i += j;
        
        if (i < num_hsps && (gapped_results[i]->context / 2) == qid)
        {
	    k = gapped_results[i]->query_gapped_start;
	    j += k;
	    i += k;
        }
		local_queries.query_offsets[qid].num_alignments = j;
    }
	
	Int4 num_queries = local_queries.GetNumQueries();
	for (i = 0; i < num_queries; ++i)
	{
		HSP** hsps = (HSP**)gapped_results.get_data() + local_queries.query_offsets[i].result_offset;
		j = local_queries.query_offsets[i].num_alignments;
		result_writter->PrintOneResult(local_queries, i, hsps, j, sbp);
	}
}

void SearchWorker::CleanUp()
{
    Int4 num_hsps = gapped_results.size();
    Int4 i;
    for (i = 0; i < num_hsps; ++i)
    {
        HSPDelete(soa, gapped_results[i]);
        gapped_results[i] = NULL;
    }
    
    if (sbp) { delete sbp; sbp = NULL; }

    if (word_params) { delete word_params; word_params = NULL; }
    if (ext_params) { delete ext_params; ext_params = NULL; }
    if (hit_params) { delete hit_params; hit_params = NULL; }
    if (eff_params)  { delete eff_params; eff_params = NULL; }
    if (score_params) { delete score_params;  score_params = NULL; }
    
    if (gapped_aligner) { delete gapped_aligner; gapped_aligner = NULL; }
    
    eff_options->db_length = 0;
    eff_options->dbseq_num = 0;
    eff_options->num_searchspaces = 1;
    eff_options->searchsp_eff = NULL;
    
    soa.Clear();
	results->clear();
}

void SetupThreadQueries(QueryInfo& qin,
                        QueryInfo& qout,
                        int qid,
                        int nthreads);

void SearchWorker::Go()
{	
	local_queries.Clear();
	SetupThreadQueries(*global_queries, local_queries, thread_id, options->running_options->num_threads);
	if (local_queries.max_length == 0) return;
	
    local_queries.MakeBlastnaQuery();
    local_queries.GetScanRanges(options->seed_options->seed_size, dust_masker, window_masker);
    BLAST_GapAlignSetup();
	
    Seeding();

	PrelimSearchStage(seeds,
					  &local_queries,
					  dbinfo,
					  options->seed_options->seed_size,
					  ungapped_results,
					  gapped_results,
					  word_params,
					  sbp,
					  hit_params,
					  &itree,
					  options->hit_options->min_diag_separation,
					  gapped_aligner,
					  soa);

    BLAST_ComputeTraceback(); 
	
    OutputResult(results);

	local_queries.query_names.set_data(NULL, 0, 0, TRUE);
    local_queries.query_offsets.set_data(NULL, 0, 0, TRUE);
    local_queries.org_queries.set_data(NULL, 0, 0, TRUE);
}

void SearchWorker::FlushResults(FILE* file)
{
	results->FlushResults(file);
}
