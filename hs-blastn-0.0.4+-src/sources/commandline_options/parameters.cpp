#include "options.h"
#include "parameters.h"
#include "stat_functions.h"

#include <cmath>

BlastExtensionParameters::BlastExtensionParameters(
                            BlastExtensionOptions* ex_options,
                            BlastScoreBlk* sbp,
                            QueryInfo* query_info)
{
    ASSERT(sbp->kbp != NULL);

    if (sbp->kbp)
    {
        Blast_KarlinBlk* kbp = NULL;
        Int2 status = 0;
        status = s_BlastFindValidKarlinBlk(sbp->kbp, query_info, &kbp);
        if (status != 0) return;
    }

    options = ex_options;

    if (sbp->kbp_gap)
    {
        double min_lambda = s_BlastFindSmallestLambda(sbp->kbp_gap, query_info, NULL);
        gap_x_dropoff = (Int4)(options->gap_x_dropoff * NCBIMATH_LN2 / min_lambda);
        gap_x_dropoff_final = (Int4)
            MAX(options->gap_x_dropoff_final * NCBIMATH_LN2 / min_lambda, gap_x_dropoff);
    }
}

BlastExtensionParameters::~BlastExtensionParameters()
{

}

BlastHitSavingParameters::BlastHitSavingParameters(
                             BlastHitSavingOptions* hs_options,
                             const BlastScoreBlk* sbp,
                             const QueryInfo* query_info,
                             Int8 avg_subj_length)
{
    Boolean gapped_calculation = TRUE;
    Int2 status = 0;

    ASSERT(hs_options != NULL);
    ASSERT(sbp != NULL);

    if (!sbp->kbp_gap)
    {
        gapped_calculation = FALSE;
    }

    mask_level = 101;
    do_sum_stats = hs_options->do_sum_stats;
    options = hs_options;
    cutoffs = (BlastGappedCutoffs*)calloc(
                                query_info->Size() * 2,
                                sizeof(BlastGappedCutoffs));

    if (options->low_score_perc > 0.00001)
    {
        low_score = (int*)calloc(query_info->Size(), sizeof(int));
    }
    else
    {
        low_score = NULL;
    }

    status = BlastHitSavingParametersUpdate(sbp, query_info, avg_subj_length);
    ASSERT(status == 0);
}

Int2
BlastHitSavingParameters::BlastHitSavingParametersUpdate(
                            const BlastScoreBlk* sbp,
                            const QueryInfo* query_info,
                            Int8 avg_subject_length)
{
    Blast_KarlinBlk** kbp_array;
    double scale_factor = sbp->scale_factor;
    Boolean gapped_calculation = TRUE;
    Int4 context;

    ASSERT(query_info);

    if (sbp->kbp_gap)
    {
        kbp_array = sbp->kbp_gap;
    }
    else if (sbp->kbp)
    {
        kbp_array = sbp->kbp;
        gapped_calculation = FALSE;
    }
    else
    {
        return -1;
    }

    if (options->mask_level >= 0)
    {
        mask_level = options->mask_level;
    }

    if (options->cutoff_score > 0)
    {
        Int4 new_cutoff = options->cutoff_score * (Int4)sbp->scale_factor;
        for (context = 0; context < query_info->Size() * 2; ++context)
        {
            cutoffs[context].cutoff_score = new_cutoff;
            cutoffs[context].cutoff_score_max = new_cutoff;
        }
        cutoff_score_min = new_cutoff;
    } else
    {
        Int4 cutoff_min = INT4_MAX;
        Blast_KarlinBlk* kbp;

        for (context = 0; context < query_info->Size() * 2; ++context)
        {
            Int8 searchsp;
            Int4 new_cutoff = 1;
            double evalue = options->expect_value;
			
			if (!(query_info->contexts[context].is_valid))
			{
				cutoffs[context].cutoff_score = INT4_MAX;
				continue;
			}

            kbp = kbp_array[context];
            if (kbp == NULL) continue;
            ASSERT(s_BlastKarlinBlkIsValid(kbp));
            searchsp = sbp->eff_searchsp[context];

            BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
            cutoffs[context].cutoff_score = new_cutoff;
            cutoffs[context].cutoff_score_max = new_cutoff;
        }

        for (context = 0; context < query_info->Size() * 2; ++context)
        {
            if (query_info->contexts[context].is_valid && sbp->kbp[context] != NULL)
            {
                cutoffs[context].cutoff_score *= (Int4)scale_factor;
                cutoffs[context].cutoff_score_max *= (Int4)scale_factor;
                cutoff_min = MIN(cutoff_min, cutoffs[context].cutoff_score);
            }
        }
        cutoff_score_min = cutoff_min;
    }

    return 0;
}

BlastHitSavingParameters::~BlastHitSavingParameters()
{
    sfree(cutoffs);
}

BlastInitialWordParameters::BlastInitialWordParameters(
		const BlastInitialWordOptions* word_options,
		const BlastHitSavingParameters* hit_params,
		const BlastScoreBlk* sbp,
	        QueryInfo* query_info,
		Int8 subject_length) : cutoffs(NULL)
{
	Int2 status;
	Blast_KarlinBlk* kbp;
	Int4 i;

	ASSERT(word_options);
	ASSERT(sbp);
	status = s_BlastFindValidKarlinBlk(sbp->kbp, query_info, &kbp);
	if (status != 0) return;

	ungapped_extension = TRUE;
	cutoffs = (BlastUngappedCutoffs*)calloc(sizeof(BlastUngappedCutoffs),
			query_info->Size() * 2);
	ASSERT(cutoffs != NULL);
	options = word_options;

	for (i = 0; i < query_info->Size() * 2; ++i)
	{
		if (!(query_info->contexts[i].is_valid)) continue;
		kbp = sbp->kbp[i];
		ASSERT(s_BlastKarlinBlkIsValid(kbp));
		cutoffs[i].x_dropoff_init = (Int4)
				(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));

	}

	status = BlastInitialWordParamtersUpdate(hit_params, sbp, query_info, subject_length);

	Int4 reward = sbp->reward;
	Int4 penalty = sbp->penalty;
	Int4* table = nucl_score_table;

	for (i = 0; i < 256; ++i)
	{
		Int4 score = 0;
        if (i & 3) score += penalty; else score += reward;
        if ((i >> 2) & 3) score += penalty; else score += reward;
        if ((i >> 4) & 3) score += penalty; else score += reward;
        if (i >> 6) score += penalty; else score += reward;
        table[i] = score;
	}
}

Int2 BlastInitialWordParameters::BlastInitialWordParamtersUpdate(
		const BlastHitSavingParameters* hit_params,
		const BlastScoreBlk* sbp,
		QueryInfo* query_info,
		Int8 subj_length)
{
	Blast_KarlinBlk** kbp_array;
	Boolean gapped_calculation = TRUE;
	double gap_decay_rate = 0.0;
	Int4 cutoff_min = INT4_MAX;
	Int4 xdrop_max = 0;
	Int4 context;

	ASSERT(sbp);
	ASSERT(hit_params);
	ASSERT(query_info);

	if (sbp->kbp_gap)
		kbp_array = sbp->kbp_gap;
	else if (sbp->kbp)
	{
		kbp_array = sbp->kbp;
		gapped_calculation = FALSE;
	}
	else return -1;

	for (context = 0; context < query_info->Size() * 2; ++context)
	{
		Int4 gap_trigger = INT4_MAX;
		Blast_KarlinBlk* kbp;
		Int4 new_cutoff = 1;
		BlastUngappedCutoffs* curr_cutoffs = cutoffs + context;
		
		if (!(query_info->contexts[context].is_valid))
		{
			curr_cutoffs->cutoff_score = INT4_MAX;
			continue;
		}

		if (sbp->kbp)
		{
			kbp = sbp->kbp[context];
			if (s_BlastKarlinBlkIsValid(kbp))
			{
				gap_trigger = (Int4)((options->gap_trigger * NCBIMATH_LN2 +
						kbp->logK) / kbp->Lambda);
			}
		}

		new_cutoff = gap_trigger;
		new_cutoff = (Int4)(sbp->scale_factor * new_cutoff);
                new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
		curr_cutoffs->cutoff_score = new_cutoff;

		if (curr_cutoffs->x_dropoff_init == 0)
			curr_cutoffs->x_dropoff = new_cutoff;
		else
			curr_cutoffs->x_dropoff = curr_cutoffs->x_dropoff_init;

		if (new_cutoff < cutoff_min)
			cutoff_min = new_cutoff;

		if (xdrop_max < curr_cutoffs->x_dropoff)
			xdrop_max = curr_cutoffs->x_dropoff;

		curr_cutoffs->reduced_nucl_cutoff_score =
				(Int4)(0.9 * new_cutoff);

		if (0)
		{
			fprintf(stderr, "[%s, line %d] context: %d, query length: %d\n",
					__func__, __LINE__, context, query_info->GetSeqLength(context));
			fprintf(stderr,
					"[%s, line %d] x_drop_init: %d, x_dropoff: %d, cutoff_score: %d\n",
					__func__, __LINE__, curr_cutoffs->x_dropoff_init,
					curr_cutoffs->x_dropoff, curr_cutoffs->cutoff_score);
		}
	}

	cutoff_score_min = cutoff_min;
	x_dropoff_max = xdrop_max;

	return 0;
}

BlastInitialWordParameters::~BlastInitialWordParameters()
{
	if (cutoffs) free(cutoffs);
}
