/*  
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Tom Madden
 *
 */

/** @file blast_stat.h
 * Definitions and prototypes used by blast_stat.c to calculate BLAST
 * statistics. @todo FIXME: needs doxygen comments
 */

#ifndef BLAST_STAT_H
#define	BLAST_STAT_H

#include "def.h"
#include "options.h"

struct Blast_KarlinBlk
{
    double Lambda;
    double K;
    double logK;
    double H;
};

/**< minimum allowed score (for one letter comparison). */
#define BLAST_SCORE_MIN INT2_MIN   
/**< maximum allowed score (for one letter comparison). */
#define BLAST_SCORE_MAX INT2_MAX   
/**< maximum allowed range of BLAST scores. */
#define BLAST_SCORE_RANGE_MAX   (BLAST_SCORE_MAX - BLAST_SCORE_MIN) 

struct Blast_ScoreFreq
{
    Int4 score_min;
    Int4 score_max;
    Int4 obs_min;
    Int4 obs_max;
    double score_avg;
    double* sprob0;
    double* sprob;
};

struct SBlastScoreMatrix
{
    int** data;
    int ncols;
    int nrows;
    double* freqs;
    double lambda;
};

/**
Stores the letter frequency of a sequence or database.
*/
struct Blast_ResFreq {
    Uint1   alphabet_code;    /**< indicates alphabet. */
    double* prob;       /**< letter probs, (possible) non-zero offset. */
    double* prob0;            /**< probs, zero offset. */
};

#include "query_info.h"

struct BlastEffectiveLengthsParameters;

struct BlastScoreBlk
{
    Uint1 alphabet_code;
    Int2  alphabet_size;
    Int2 alphabet_start;
    SBlastScoreMatrix* matrix;
    Int4 loscore;
    Int4 hiscore;
    Int4 penalty;
    Int4 reward;
    double scale_factor;
    Blast_ScoreFreq** sfp;
    Blast_KarlinBlk** kbp;
    Blast_KarlinBlk** kbp_gap;
    Blast_KarlinBlk*  kbp_ideal;
    Int4 number_of_contexts;
    Uint1* ambiguous_res;
    Int2 ambig_size, ambig_occupy;
    Boolean round_down;

    QueryInfo* query_blk;
    Int8* eff_searchsp;
    Int4* length_adjustment;

    BlastScoreBlk(Uint1 alphabet, QueryInfo* queries);
    Int2 ScoreBlkInit(const BlastScoringOptions* scoring_options, double sf);
    Int2 ScoreBlkKbpUngappedCalc();
    Int2 ScoreBlkKbpGappedCalc(const BlastScoringOptions* scoring_options);
    Int2 ScoreBlkKbpIdealCalc();
    Int2 ScoreBlkMatrixInit(const BlastScoringOptions* scoring_options);
    Int2 ScoreBlkMatrixFill();
    Int2 ScoreBlkSetAmbigRes(char ambig_res);

    Int2 BLAST_CalcEffLength(const BlastScoringOptions* scoring_options,
    		const BlastEffectiveLengthsParameters* eff_len_params,
    		QueryInfo* query_info);
    Int8 s_GetEffectiveSearchSpaceForContext(
    		const BlastEffectiveLengthsOptions* eff_len_options,
    		int context_index);

    void ScoreBlkFree();
    ~BlastScoreBlk();
};

class UngappedBlastnKarlinBlock
{
public:
    Int2 KarlinBlkUngappedCalc(Blast_KarlinBlk* kbp, Blast_ScoreFreq* sfp);

private:
    double KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess);
    double NlmKarlinLambdaNR(double* probs, Int4 d, Int4 low, Int4 high, double lambda0,
		      double tolx, Int4 itmax, Int4 maxNewton, Int4* itn);
    double KarlinLtoH(Blast_ScoreFreq* sfp, double lambda);
    double KarlinLHtoK(Blast_ScoreFreq* sfp, double lambda, double H);
};

/** Number of statistical parameters in each row of the precomputed tables. */
/**< originally 8, now 11 to support Spouge's FSC. see notes below */
#define BLAST_NUM_STAT_VALUES 11  

/** Holds values (gap-opening, extension, etc.) for a matrix. */
typedef double array_of_8[BLAST_NUM_STAT_VALUES];

class GappedBlastnKarlinBlock
{
public:
    Int2 KarlinBlkNuclGappedCalc(Blast_KarlinBlk* kbp, Int4 gap_open,
            Int4 gap_extend, Int4 reward, Int4 penalty,
            Blast_KarlinBlk* kbp_ungap, Boolean* round_down);

    Int2 SplitArrayOf8(const array_of_8* input, const array_of_8** normal,
            const array_of_8** non_affine, Boolean* split);
    Int2 GetNuclValuesArray(Int4 reward, Int4 penalty, Int4* array_size,
            array_of_8** normal, array_of_8** non_affine,
            Int4* gap_open_max, Int4* gap_extend_max, Boolean* round_down);
    Int2 AdjustGapParametersByGcd(array_of_8* normal, array_of_8* linear,
            int size, Int4* gap_existence_max, Int4* gap_extend_max, int divisor);
};

#endif	/* BLAST_STAT_H */

