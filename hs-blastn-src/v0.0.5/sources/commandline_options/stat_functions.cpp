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
 * Author: Tom Madden
 *
 */

/** @file blast_stat.c
 * Functions to calculate BLAST probabilities etc.
 * Detailed Contents:
 *
 * - allocate and deallocate structures used by BLAST to calculate
 * probabilities etc.
 *
 * - calculate residue frequencies for query and "average" database.
 *
 * - read in matrix or load it from memory.
 *
 *  - calculate sum-p from a collection of HSP's, for both the case
 *   of a "small" gap and a "large" gap, when give a total score and the
 *   number of HSP's.
 *
 * - calculate expect values for p-values.
 *
 * - calculate pseuod-scores from p-values.
 *
 */

#include "stat_functions.h"
#include "parameters.h"

#include <cstdlib>
#include <cstring>
#include <cmath>

static void**
_PSIDeallocateMatrix(void** matrix, unsigned int ncols)
{
    unsigned int i = 0;

    if (!matrix)
        return NULL;

    for (i = 0; i < ncols; i++) {
        sfree(matrix[i]);
    }
    free(matrix);
    matrix = NULL;
    return NULL;
}

SBlastScoreMatrix*
SBlastScoreMatrixFree(SBlastScoreMatrix* matrix)
{
    if ( !matrix ) {
        return NULL;
    }

    if (matrix->data) {
        matrix->data = (int**) _PSIDeallocateMatrix((void**) matrix->data,
                                                    matrix->ncols);
    }

    /* Deallocate the matrix frequencies which is used by the
     * nucleotide custom matrix reader. -RMH-
     */
    if ( matrix->freqs )
      sfree(matrix->freqs);

    free(matrix);
    matrix = NULL;
    return NULL;
}

void**
_PSIAllocateMatrix(unsigned int ncols, unsigned int nrows,
                   unsigned int data_type_sz)
{
    void** retval = NULL;
    unsigned int i = 0;

    retval = (void**) malloc(sizeof(void*) * ncols);
    if ( !retval ) {
        return NULL;
    }

    for (i = 0; i < ncols; i++) {
        retval[i] = (void*) calloc(nrows, data_type_sz);
        if ( !retval[i] ) {
            retval = _PSIDeallocateMatrix(retval, i);
            break;
        }
    }
    return retval;
}

SBlastScoreMatrix*
SBlastScoreMatrixNew(size_t ncols, size_t nrows)
{
    SBlastScoreMatrix* retval = NULL;

    retval = (SBlastScoreMatrix*) calloc(1, sizeof(SBlastScoreMatrix));
    if ( !retval ) {
        return SBlastScoreMatrixFree(retval);
    }

    retval->data = (int**) _PSIAllocateMatrix(ncols, nrows, sizeof(int));
    if ( !retval->data ) {
        return SBlastScoreMatrixFree(retval);
    }

    /* Allocate additional attributes for use with custom
     * nucleotide matrices. -RMH-
     */
    retval->freqs = (double *) calloc(ncols, sizeof(double));
    retval->lambda = 0;

    retval->ncols = ncols;
    retval->nrows = nrows;
    return retval;
}

Blast_ScoreFreq*
Blast_ScoreFreqFree(Blast_ScoreFreq* sfp)
{
   if (sfp == NULL)
      return NULL;

   if (sfp->sprob0 != NULL)
      sfree(sfp->sprob0);
   free(sfp);
   sfp = NULL;
   return sfp;
}

Blast_KarlinBlk*
Blast_KarlinBlkFree(Blast_KarlinBlk* kbp)

{
   free(kbp);
   kbp = NULL;

   return kbp;
}

const Uint1 IUPACNA_TO_BLASTNA[128]={
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};

const Uint1 IUPACNA_TO_NCBI4NA[128]={
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
 0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const Uint1 NCBI4NA_TO_BLASTNA[BLASTNA_SIZE] = {
    15,     /* Gap, 0 */
     0,     /* A,   1 */
     1,     /* C,   2 */
     6,     /* M,   3 */
     2,     /* G,   4 */
     4,     /* R,   5 */
     9,     /* S,   6 */
    13,     /* V,   7 */
     3,     /* T,   8 */
     8,     /* W,   9 */
     5,     /* Y,  10 */
    12,     /* H,  11 */
     7,     /* K,  12 */
    11,     /* D,  13 */
    10,     /* B,  14 */
    14      /* N,  15 */
};

const Uint1 BLASTNA_TO_NCBI4NA[BLASTNA_SIZE] = {
     1,     /* A, 0 */
     2,     /* C, 1 */
     4,     /* G, 2 */
     8,     /* T, 3 */
     5,     /* R, 4 */
    10,     /* Y, 5 */
     3,     /* M, 6 */
    12,     /* K, 7 */
     9,     /* W, 8 */
     6,     /* S, 9 */
    14,     /* B, 10 */
    13,     /* D, 11 */
    11,     /* H, 12 */
     7,     /* V, 13 */
    15,     /* N, 14 */
     0      /* Gap, 15 */
};

long BLAST_Nint(double x)
{
   x += (x >= 0. ? 0.5 : -0.5);
   return (long)x;
}

/** Fill in the matrix for blastn using the penaly and rewards
 * The query sequence alphabet is blastna, the subject sequence
 * is ncbi2na.  The alphabet blastna is defined in blast_stat.h
 * and the first four elements of blastna are identical to ncbi2na.
 * @param sbp the BlastScoreBlk on which reward, penalty, and matrix will be set [in|out]
 * @return zero on success.
*/

Int2 BlastScoreBlkNuclMatrixCreate(BlastScoreBlk* sbp)
{
    Int2  index1, index2, degen;
    Int2 degeneracy[BLASTNA_SIZE+1];
    Int4 reward; /* reward for match of bases. */
    Int4 penalty; /* cost for mismatch of bases. */
    Int4** matrix; /* matrix to be populated. */
    /* How many of the first bases are ambiguous (four, of course). */
    const int k_number_non_ambig_bp = 4;

    ASSERT(sbp);
    ASSERT(sbp->alphabet_size == BLASTNA_SIZE);
    ASSERT(sbp->matrix);
    ASSERT(sbp->matrix->ncols == BLASTNA_SIZE);
    ASSERT(sbp->matrix->nrows == BLASTNA_SIZE);

    reward = sbp->reward;
    penalty = sbp->penalty;
    matrix = sbp->matrix->data;

    for (index1 = 0; index1<BLASTNA_SIZE; index1++)
        for (index2 = 0; index2<BLASTNA_SIZE; index2++)
            matrix[index1][index2] = 0;

    /* In blastna the 1st four bases are A, C, G, and T, exactly as it is
       ncbi2na. */
    /* ncbi4na gives them the value 1, 2, 4, and 8.  */
    /* Set the first four bases to degen. one */
    for (index1=0; index1<k_number_non_ambig_bp; index1++)
        degeneracy[index1] = 1;

    for (index1=k_number_non_ambig_bp; index1<BLASTNA_SIZE; index1++) {
        degen=0;
        for (index2=0; index2<k_number_non_ambig_bp; index2++) /* ncbi2na */
        {
            if (BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2])
                degen++;
        }
        degeneracy[index1] = degen;
    }


    for (index1=0; index1<BLASTNA_SIZE; index1++) {
        for (index2=index1; index2<BLASTNA_SIZE; index2++) {
            if (BLASTNA_TO_NCBI4NA[index1] & BLASTNA_TO_NCBI4NA[index2]) {
                /* round up for positive scores, down for negatives. */
                matrix[index1][index2] =
                    BLAST_Nint( (double) ((degeneracy[index2]-1)*penalty +
                                          reward)/ (double) degeneracy[index2]);
                if (index1 != index2)
                {
                      matrix[index2][index1] = matrix[index1][index2];
                }
            }
            else
            {
                matrix[index1][index2] = penalty;
                matrix[index2][index1] = penalty;
            }
        }
    }

    /* The value of 15 is a gap, which is a sentinel between strands in
    the ungapped extension algorithm */
    for (index1=0; index1<BLASTNA_SIZE; index1++)
        matrix[BLASTNA_SIZE-1][index1] = INT4_MIN / 2;
    for (index1=0; index1<BLASTNA_SIZE; index1++)
        matrix[index1][BLASTNA_SIZE-1] = INT4_MIN / 2;

    return 0;
}

/** Sets maximum and minimum scores on the BlastScoreBlk for a
 * given matrix
 * @param sbp the BlastScoreBlk on which loscore and hiscore
 *   will be set [in|out]
 * @return zero on success
 */
Int2
BlastScoreBlkMaxScoreSet(BlastScoreBlk* sbp)
{
    Int4 score;
    Int4 ** matrix;
    Int2 index1, index2;

    sbp->loscore = BLAST_SCORE_MAX;
    sbp->hiscore = BLAST_SCORE_MIN;
    matrix = sbp->matrix->data;
    for (index1=0; index1<sbp->alphabet_size; index1++)
    {
      for (index2=0; index2<sbp->alphabet_size; index2++)
      {
         score = matrix[index1][index2];
         if (score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX)
            continue;
         if (sbp->loscore > score)
            sbp->loscore = score;
         if (sbp->hiscore < score)
            sbp->hiscore = score;
      }
    }
    /* If the lo/hi-scores are BLAST_SCORE_MIN/BLAST_SCORE_MAX, (i.e., for
    gaps), then use other scores. */

    if (sbp->loscore < BLAST_SCORE_MIN)
      sbp->loscore = BLAST_SCORE_MIN;
    if (sbp->hiscore > BLAST_SCORE_MAX)
      sbp->hiscore = BLAST_SCORE_MAX;

    return 0;
}

Int2
Blast_ScoreBlkMatrixFill(BlastScoreBlk* sbp)
{
    Boolean matrix_found = FALSE;
    Int2 status = 0;

    /* For nucleotide case we first create a default matrix, based on the
       match and mismatch scores. */
    if (sbp->alphabet_code == BLASTNA_SEQ_CODE) {
		if ((status = BlastScoreBlkNuclMatrixCreate(sbp)) != 0)
			return status;
		matrix_found = TRUE;
    }

    if (matrix_found == FALSE)
        return -1;

    if ( (status=BlastScoreBlkMaxScoreSet(sbp)) != 0)
         return status;

    return status;
}

Blast_ResFreq*
Blast_ResFreqFree(Blast_ResFreq* rfp)
{
   if (rfp == NULL)
      return NULL;

   if (rfp->prob0 != NULL)
      sfree(rfp->prob0);

   free(rfp);
   rfp = NULL;

   return rfp;
}

Blast_ResFreq*
Blast_ResFreqNew(const BlastScoreBlk* sbp)
{
   Blast_ResFreq* rfp;

   if (sbp == NULL)
   {
      return NULL;
   }

   rfp = (Blast_ResFreq*) calloc(1, sizeof(Blast_ResFreq));
   if (rfp == NULL)
      return NULL;

   rfp->alphabet_code = sbp->alphabet_code;

   rfp->prob0 = (double*) calloc(sbp->alphabet_size, sizeof(double));
   if (rfp->prob0 == NULL)
   {
      rfp = Blast_ResFreqFree(rfp);
      return rfp;
   }
   rfp->prob = rfp->prob0 - sbp->alphabet_start;

   return rfp;
}

struct BLAST_LetterProb
{
	char ch;
	double p;
};

static BLAST_LetterProb nt_prob[] = {
		{ 'A', 25.00 },
		{ 'C', 25.00 },
		{ 'G', 25.00 },
		{ 'T', 25.00 }
};

/** Normalizes all the residue frequencies and then normalizes them to "norm".
 * If "norm" is one, then they will all sum to one.
 * @param sbp needed for alphabet information [in]
 * @param rfp array of residue frequencies to be normalized [in|out]
 * @param norm value to normalize to [in]
 * @return zero on success, 1 otherwise
*/
static Int2
Blast_ResFreqNormalize(const BlastScoreBlk* sbp, Blast_ResFreq* rfp, double norm)
{
   Int2  alphabet_stop, index;
   double   sum = 0., p;

   if (norm == 0.)
      return 1;

   alphabet_stop = sbp->alphabet_start + sbp->alphabet_size;
   for (index=sbp->alphabet_start; index<alphabet_stop; index++)
   {
      p = rfp->prob[index];
      if (p < 0.)
         return 1;
      sum += p;
   }
   if (sum <= 0.)
      return 0;

   for (index=sbp->alphabet_start; index<alphabet_stop; index++)
   {
      rfp->prob[index] /= sum;
      rfp->prob[index] *= norm;
   }
   return 0;
}

Int2
Blast_ResFreqStdComp(const BlastScoreBlk* sbp, Blast_ResFreq* rfp)
{
        Uint4 index;

   {  /* beginning of blastna and ncbi2na are the same. */
      /* Only blastna used  for nucleotides. */
      for (index=0; index<DIM(nt_prob); index++)
      {
         rfp->prob[index] = nt_prob[index].p;
      }
   }

   Blast_ResFreqNormalize(sbp, rfp, 1.0);

   return 0;
}

struct Blast_ResComp
{
	Uint1 alphabet_code;
	Int4* comp;
	Int4* comp0;
};

static Blast_ResComp*
BlastResCompDestruct(Blast_ResComp* rcp)
{
   if (rcp == NULL)
      return NULL;

   if (rcp->comp0 != NULL)
      sfree(rcp->comp0);

   free(rcp);
   rcp = NULL;
   return NULL;
}

static Blast_ResComp*
BlastResCompNew(const BlastScoreBlk* sbp)
{
   Blast_ResComp* rcp;

   rcp = (Blast_ResComp*) calloc(1, sizeof(Blast_ResComp));
   if (rcp == NULL)
      return NULL;

   rcp->alphabet_code = sbp->alphabet_code;

/* comp0 has zero offset, comp starts at 0, only one
array is allocated.  */
   rcp->comp0 = (Int4*) calloc(sbp->alphabet_size, sizeof(Int4));
   if (rcp->comp0 == NULL)
   {
      rcp = BlastResCompDestruct(rcp);
      return rcp;
   }

   rcp->comp = rcp->comp0 - sbp->alphabet_start;

   return rcp;
}

/** Store the composition of a (query) string.
 * @param sbp needed for alphabet information [in]
 * @param rcp object to be filled in [in|out]
 * @param str sequence to have composition calculated for [in]
 * @param length length of sequence [in]
 * @return zero on success, 1 otherwise.
*/
static Int2
BlastResCompStr(const BlastScoreBlk* sbp, Blast_ResComp* rcp, char* str, Int4 length)
{
   char* lp,* lpmax;
   Int2 index;
        Uint1 mask;

   if (sbp == NULL || rcp == NULL || str == NULL)
      return 1;

   if (rcp->alphabet_code != sbp->alphabet_code)
      return 1;

        /* For megablast, check only the first 4 bits of the sequence values */
        mask = (0x0f);

/* comp0 starts at zero and extends for "num", comp is the same array, but
"start_at" offset. */
   for (index=0; index<(sbp->alphabet_size); index++)
      rcp->comp0[index] = 0;

   for (lp = str, lpmax = lp+length; lp < lpmax; lp++)
   {
      ++rcp->comp[(int)(*lp & mask)];
   }

   /* Don't count ambig. residues. */
   for (index=0; index<sbp->ambig_occupy; index++)
   {
      rcp->comp[sbp->ambiguous_res[index]] = 0;
   }

   return 0;
}

/** Sets prob elements of Blast_ResFreq to zero
 * @param sbp needed for alphabet information [in]
 * @param rfp contains elements to be zeroed [in|out]
 * @return zero on success.
 */
static Int2
Blast_ResFreqClr(const BlastScoreBlk* sbp, Blast_ResFreq* rfp)
{
   Int2  alphabet_max, index;

   alphabet_max = sbp->alphabet_start + sbp->alphabet_size;
   for (index=sbp->alphabet_start; index<alphabet_max; index++)
                rfp->prob[index] = 0.0;

        return 0;
}

/** Calculate the residue frequencies associated with the provided ResComp
 *  This function takes into account the composition of a given sequence
 *  (expressed through rcp) rather than just doing it for a standard distribution.
 * @param sbp contains alphabet information [in]
 * @param rfp object to be filled in [in|out]
 * @param rcp object with composition information [in]
 * @return zero on success, 1 on failure
*/
static Int2
Blast_ResFreqResComp(const BlastScoreBlk* sbp, Blast_ResFreq* rfp,
                     const Blast_ResComp* rcp)
{
   Int2  alphabet_max, index;
   double   sum = 0.;

   if (rfp == NULL || rcp == NULL)
      return 1;

   if (rfp->alphabet_code != rcp->alphabet_code)
      return 1;

   alphabet_max = sbp->alphabet_start + sbp->alphabet_size;
   for (index=sbp->alphabet_start; index<alphabet_max; index++)
      sum += rcp->comp[index];

   if (sum == 0.) {
      Blast_ResFreqClr(sbp, rfp);
      return 0;
   }

   for (index=sbp->alphabet_start; index<alphabet_max; index++)
      rfp->prob[index] = rcp->comp[index] / sum;

   return 0;
}

/** Fills in residue frequences for a given sequence.
 * @param sbp needed for alphabet information [in]
 * @param rfp object to be populated [in|out]
 * @param string sequence for calculation [in]
 * @param length length of above sequence [in]
 */
Int2
Blast_ResFreqString(const BlastScoreBlk* sbp, Blast_ResFreq* rfp, char* string, Int4 length)
{
   Blast_ResComp* rcp;

   rcp = BlastResCompNew(sbp);

   BlastResCompStr(sbp, rcp, string, length);

   Blast_ResFreqResComp(sbp, rfp, rcp);

   rcp = BlastResCompDestruct(rcp);

   return 0;
}

/** Check that the lo and hi score are within the allowed ranges
 * @param lo the lowest permitted value [in]
 * @param hi the highest permitted value [in]
 * @return zero on success, 1 otherwise
 */

static Int2
BlastScoreChk(Int4 lo, Int4 hi)
{
   if (lo >= 0 || hi <= 0 ||
         lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX)
      return 1;

   if (hi - lo > BLAST_SCORE_RANGE_MAX)
      return 1;

   return 0;
}

Blast_ScoreFreq*
Blast_ScoreFreqNew(Int4 score_min, Int4 score_max)
{
   Blast_ScoreFreq*  sfp;
   Int4  range;

   if (BlastScoreChk(score_min, score_max) != 0)
      return NULL;

   sfp = (Blast_ScoreFreq*) calloc(1, sizeof(Blast_ScoreFreq));
   if (sfp == NULL)
      return NULL;

   range = score_max - score_min + 1;
   sfp->sprob = (double*) calloc(range, sizeof(double));
   if (sfp->sprob == NULL)
   {
      Blast_ScoreFreqFree(sfp);
      return NULL;
   }

   sfp->sprob0 = sfp->sprob;
   sfp->sprob -= score_min;        /* center around 0 */
   sfp->score_min = score_min;
   sfp->score_max = score_max;
   sfp->obs_min = sfp->obs_max = 0;
   sfp->score_avg = 0.0;
   return sfp;
}

/** Calculates the score frequencies.
 *
 * @param sbp object with scoring information [in]
 * @param sfp object to hold frequency information [in|out]
 * @param rfp1 letter frequencies for first sequence (query) [in]
 * @param rfp2 letter frequencies for second sequence (database) [in]
 * @return zero on success
 */
Int2
BlastScoreFreqCalc(const BlastScoreBlk* sbp, Blast_ScoreFreq* sfp, Blast_ResFreq* rfp1, Blast_ResFreq* rfp2)
{
   Int4 **  matrix;
   Int4  score, obs_min, obs_max;
   double      score_sum, score_avg;
   Int2     alphabet_start, alphabet_end, index1, index2;

   if (sbp == NULL || sfp == NULL)
      return 1;

   if (sbp->loscore < sfp->score_min || sbp->hiscore > sfp->score_max)
      return 1;

   for (score = sfp->score_min; score <= sfp->score_max; score++)
      sfp->sprob[score] = 0.0;

   matrix = sbp->matrix->data;

   alphabet_start = sbp->alphabet_start;
   alphabet_end = alphabet_start + sbp->alphabet_size;
   for (index1=alphabet_start; index1<alphabet_end; index1++)
   {
      for (index2=alphabet_start; index2<alphabet_end; index2++)
      {
         score = matrix[index1][index2];
         if (score >= sbp->loscore)
         {
            sfp->sprob[score] += rfp1->prob[index1] * rfp2->prob[index2];
         }
      }
   }

   score_sum = 0.;
   obs_min = obs_max = BLAST_SCORE_MIN;
   for (score = sfp->score_min; score <= sfp->score_max; score++)
   {
      if (sfp->sprob[score] > 0.)
      {
         score_sum += sfp->sprob[score];
         obs_max = score;
         if (obs_min == BLAST_SCORE_MIN)
            obs_min = score;
      }
   }
   sfp->obs_min = obs_min;
   sfp->obs_max = obs_max;

   score_avg = 0.0;
   if (score_sum > 0.0001 || score_sum < -0.0001)
   {
      for (score = obs_min; score <= obs_max; score++)
      {
         sfp->sprob[score] /= score_sum;
         score_avg += score * sfp->sprob[score];
      }
   }
   sfp->score_avg = score_avg;

   return 0;
}

Blast_KarlinBlk*
Blast_KarlinBlkNew(void)

{
   Blast_KarlinBlk* kbp;

   kbp = (Blast_KarlinBlk*) calloc(1, sizeof(Blast_KarlinBlk));

   return kbp;
}

Int2 Blast_KarlinBlkCopy(Blast_KarlinBlk* kbp_to, Blast_KarlinBlk* kbp_from)
{
   if (!kbp_to || !kbp_from)
      return -1;

   kbp_to->Lambda = kbp_from->Lambda;
   kbp_to->K = kbp_from->K;
   kbp_to->logK = kbp_from->logK;
   kbp_to->H = kbp_from->H;
   return 0;
}

Int2
s_BlastFindValidKarlinBlk(Blast_KarlinBlk** kbp_in,
		const QueryInfo* query_info,
		Blast_KarlinBlk** kbp_ret)
{
	Int4 i;
	Int2 status = 1;
	ASSERT(kbp_in && query_info && kbp_ret);

	for (i = 0; i < query_info->Size() * 2; ++i)
	{
		ASSERT(s_BlastKarlinBlkIsValid(kbp_in[i]) == query_info->contexts[i].is_valid);
		if (s_BlastKarlinBlkIsValid(kbp_in[i]))
		{
			*kbp_ret = kbp_in[i];
			status = 0;
			break;
		}
	}
	return status;
}

Boolean s_BlastKarlinBlkIsValid(const Blast_KarlinBlk* kbp)
{
    if ( !kbp ) {
        return FALSE;
    } else {
        return (kbp->Lambda > 0 && kbp->K > 0 && kbp->H > 0);
    }
}

double
s_BlastFindSmallestLambda(Blast_KarlinBlk** kbp_in,
                          const QueryInfo* query_info,
                          Blast_KarlinBlk** kbp_out)
{
    Int4 i;
    double min_lambda = (double) INT4_MAX;

    ASSERT(kbp_in && query_info);

    for (i = 0; i < query_info->Size() * 2; i++) {
		ASSERT(s_BlastKarlinBlkIsValid(kbp_in[i]) == query_info->contexts[i].is_valid);
        if (s_BlastKarlinBlkIsValid(kbp_in[i])) {
            if (min_lambda > kbp_in[i]->Lambda)
            {
                min_lambda = kbp_in[i]->Lambda;
                if (kbp_out)
                  *kbp_out = kbp_in[i];
            }
        }
    }

    ASSERT(min_lambda > 0.0);
    return min_lambda;
}

static double
Powi(double x, Int4 n)
{
    double y;
    if (n == 0) return 1;
    
    if (x == 0.)
    {
        if (n < 0) return HUGE_VAL;
        return 0.;
    }
    
    if (n < 0)
    {
        x = 1./x;
        n = -n;
    }
    
    y = 1.0;
    while (n > 0)
    {
        if (n & 1) y *= x;
        n /= 2;
        x *= x;
    }
    return y;
}

/* Compute a divisor used to weight the evalue of a collection of
 * "nsegs" distinct alignments.  These divisors are used to compensate
 * for the effect of choosing the best among multiple collections of
 * alignments.  See
 *
 * Stephen F. Altschul. Evaluating the statitical significance of
 * multiple distinct local alignments. In Suhai, editior, Theoretical
 * and Computational Methods in Genome Research, pages 1-14. Plenum
 * Press, New York, 1997.
 *
 * The "decayrate" parameter of this routine is a value in the
 * interval (0,1). Typical values of decayrate are .1 and .5. */

double
BLAST_GapDecayDivisor(double decayrate, unsigned nsegs )
{
    return (1. - decayrate) * Powi(decayrate, nsegs - 1);
}

/** Calculates score from expect value and search space.
 * @param E expect value [in]
 * @param kbp contains Karlin-Altschul parameters [in]
 * @param searchsp query times database size [in]
 * @return score
 */
static Int4
BlastKarlinEtoS_simple(double E, /* Expect value */
   const Blast_KarlinBlk*  kbp,
   Int8  searchsp)   /* size of search space */
{

   double   Lambda, K, H; /* parameters for Karlin statistics */
   Int4  S;
/* Smallest float that might not cause a floating point exception in
   S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda )); below.  */
   const double kSmallFloat = 1.0e-297;

   Lambda = kbp->Lambda;
   K = kbp->K;
   H = kbp->H;
   if (Lambda < 0. || K < 0. || H < 0.0)
   {
      return BLAST_SCORE_MIN;
   }

   E = MAX(E, kSmallFloat);

   S = (Int4) (ceil( log((double)(K * searchsp / E)) / Lambda ));
   return S;
}

/*
   BlastKarlinStoE() -- given a score, return the associated Expect value

   Error return value is -1.
*/
double
BLAST_KarlinStoE_simple(Int4 S,
      Blast_KarlinBlk* kbp,
      Int8  searchsp)   /* size of search space. */
{
   double   Lambda, K, H; /* parameters for Karlin statistics */

   Lambda = kbp->Lambda;
   K = kbp->K;
   H = kbp->H;
   if (Lambda < 0. || K < 0. || H < 0.) {
      return -1.;
   }

   return (double) searchsp * exp((double)(-Lambda * S) + kbp->logK);
}

/*
   BlastCutoffs
   Calculate the cutoff score, S, and the highest expected score.
*/
Int2
BLAST_Cutoffs(Int4 *S, /* cutoff score */
   double* E, /* expected no. of HSPs scoring at or above S */
   Blast_KarlinBlk* kbp,
   Int8 searchsp, /* size of search space. */
   Boolean dodecay,  /* TRUE ==> use gapdecay feature */
   double gap_decay_rate)
{
   Int4  s = *S, es;
   double   e = *E, esave;
   Boolean     s_changed = FALSE;

   if (kbp->Lambda == -1. || kbp->K == -1. || kbp->H == -1.)
      return 1;

   /*
   Calculate a cutoff score, S, from the Expected
   (or desired) number of reported HSPs, E.
   */
   es = 1;
   esave = e;
   if (e > 0.)
   {
        if (dodecay) {
            /* Invert the adjustment to the e-value that will be applied
             * to compensate for the effect of choosing the best among
             * multiple alignments */
            if( gap_decay_rate > 0 && gap_decay_rate < 1 ) {
                e *= BLAST_GapDecayDivisor(gap_decay_rate, 1);
            }
        }
        es = BlastKarlinEtoS_simple(e, kbp, searchsp);
   }
   /*
   Pick the larger cutoff score between the user's choice
   and that calculated from the value of E.
   */
   if (es > s) {
      s_changed = TRUE;
      *S = s = es;
   }

   /*
      Re-calculate E from the cutoff score, if E going in was too high
   */
   if (esave <= 0. || !s_changed)
   {
      e = BLAST_KarlinStoE_simple(s, kbp, searchsp);
      if (dodecay) {
            /* Weight the e-value to compensate for the effect of
             * choosing the best of more than one collection of
             * distinct alignments */
            if( gap_decay_rate > 0 && gap_decay_rate < 1 ) {
                e /= BLAST_GapDecayDivisor(gap_decay_rate, 1);
            }
        }
      *E = e;
   }

   return 0;
}

Int2 BLAST_GetNucleotideGapExistenceExtendParams(Int4 reward,
                                       Int4 penalty,
                                       Int4* gap_existence,
                                       Int4* gap_extension)
{
   int array_size = 0; /* dummy parameter. */
   array_of_8* normal=NULL; /* dummy parameter */
   array_of_8* non_affine=NULL; /* dummy parameter */
   Boolean round_down = FALSE;
   int gap_existence_max=0;
   int gap_extension_max=0;
   GappedBlastnKarlinBlock gapped_karlinblk;
   Int2 status = gapped_karlinblk.GetNuclValuesArray(reward, penalty, &array_size,
		   	   	   	   	   	   	   	   	   	   	     &normal, &non_affine,
		   	   	   	   	   	   	   	   	   	   	     &gap_existence_max,
		   	   	   	   	   	   	   	   	   	   	     &gap_extension_max,
		   	   	   	   	   	   	   	   	   	   	     &round_down);

   if (status)
   {
       free(normal); normal = NULL;
       free(non_affine); non_affine = NULL;
       return status;
   }

   if (*gap_existence == 0 && *gap_extension == 0 && non_affine)
       status = 0;   /* these values are supported. */
   else
   {
         int index=0;
         Boolean found=FALSE;
         while (index < array_size)
         {
               if (*gap_existence == normal[index][0] && *gap_extension == normal[index][1])
               {
                      found = TRUE;
                      break; /* these values are supported. */
               }
               index++;
         }

         if (!found)
         {   /* If values are above max, then use. Otherwise set max values. */
             if (*gap_existence < gap_existence_max || *gap_extension < gap_extension_max)
             {
                 *gap_existence = gap_existence_max;
                 *gap_extension = gap_extension_max;
             }
         }
         status = 0;
   }
   free(normal); normal = NULL;
   free(non_affine); non_affine = NULL;
   return status;
}

Boolean BLAST_CheckRewardPenaltyScores(Int4 reward, Int4 penalty)
{
    int array_size = 0; /* dummy parameter. */
    array_of_8* normal = NULL; /* dummy parameter */
    array_of_8* non_affine = NULL; /* dummy parameter */
    Boolean round_down = FALSE;
    int gap_existence_max = 0;
    int gap_extension_max = 0;
    GappedBlastnKarlinBlock gapped_karlinblk;
    Int2 status =
    	gapped_karlinblk.GetNuclValuesArray(reward, penalty, &array_size, &normal,
                                       &non_affine, &gap_existence_max,
                                       &gap_extension_max, &round_down);

    free(normal); normal = NULL;
    free(non_affine); non_affine = NULL;

    return status == 0;
}

/** Returns the beta statistical parameter value, given the nucleotide
 * substitution scores.
 * @param reward Match reward score [in]
 * @param penalty Mismatch penalty score [in]
 * @return The value of the beta parameter.
 */
static double s_GetUngappedBeta(Int4 reward, Int4 penalty)
{
    double beta = 0;
    if ((reward == 1 && penalty == -1) ||
        (reward == 2 && penalty == -3))
        beta = -2;

    return beta;
}

Int2 Blast_GetNuclAlphaBeta(Int4 reward, Int4 penalty, Int4 gap_open,
                            Int4 gap_extend, Blast_KarlinBlk* kbp,
                            Boolean gapped_calculation,
                            double *alpha, double *beta)
{
    const int kGapOpenIndex = 0;
    const int kGapExtIndex = 1;
    const int kAlphaIndex = 5;
    const int kBetaIndex = 6;
    Int4 num_combinations = 0;
    Int4 gap_open_max = 0, gap_extend_max = 0;
    Int4 index = 0;
    array_of_8* normal=NULL;
    array_of_8* linear=NULL;
    Boolean round_down = FALSE;
    Boolean found = FALSE;
    GappedBlastnKarlinBlock gapped_karlinblk;
    Int2 status =
    		gapped_karlinblk.GetNuclValuesArray(reward,
                                       penalty,
                                       &num_combinations,
                                       &normal,
                                       &linear,
                                       &gap_open_max,
                                       &gap_extend_max,
                                       &round_down);

    if (status)
       return status;

    ASSERT(alpha && beta && kbp);

    /* For ungapped search return ungapped values of alpha and beta. */
    if (gapped_calculation && normal) {
        if (gap_open == 0 && gap_extend == 0 && linear)
        {
            *alpha = linear[0][kAlphaIndex];
            *beta = linear[0][kBetaIndex];
            found = TRUE;
        }
        else
        {

            for (index = 0; index < num_combinations; ++index) {
                if (normal[index][kGapOpenIndex] == gap_open &&
                    normal[index][kGapExtIndex] == gap_extend) {
                    *alpha = normal[index][kAlphaIndex];
                    *beta = normal[index][kBetaIndex];
                    found = TRUE;
                    break;
                }
            }
        }

    }

    /* If input values not found in tables, or if this is an ungapped search,
       return the ungapped values of alpha and beta. */
    if (!found)
    {
        *alpha = kbp->Lambda/kbp->H;
        *beta = s_GetUngappedBeta(reward, penalty);
    }

    free(linear); linear = NULL;
    free(normal); normal = NULL;
    return 0;
}

/* UngappedBlastKarlinBlock functions */

/**< K_SUMLIMIT_DEFAULT == sumlimit used in BlastKarlinLHtoK() */
#define BLAST_KARLIN_K_SUMLIMIT_DEFAULT 0.0001 

/**< LAMBDA_ACCURACY_DEFAULT == accuracy to which Lambda should be calc'd */
#define BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT    (1.e-5) 

/**< LAMBDA_ITER_DEFAULT == no. of iterations in LambdaBis = ln(accuracy)/ln(2)*/
#define BLAST_KARLIN_LAMBDA_ITER_DEFAULT        17 

/**< Initial guess for the value of Lambda in BlastKarlinLambdaNR */
#define BLAST_KARLIN_LAMBDA0_DEFAULT    0.5 

/**< upper limit on iterations for BlastKarlinLHtoK */
#define BLAST_KARLIN_K_ITER_MAX 100 

double UngappedBlastnKarlinBlock::NlmKarlinLambdaNR(double* probs, Int4 d,
						    Int4 low, Int4 high,
						    double lambda0, double tolx,
						    Int4 itmax, Int4 maxNewton, Int4* itn)
{
  Int4 k;
  double x0, x, a = 0, b = 1;
  double f = 4;  /* Larger than any possible value of the poly in [0,1] */
  Int4 isNewton = 0; /* we haven't yet taken a Newton step. */

  assert( d > 0 );

   x0 = exp( -lambda0 );
  x = ( 0 < x0 && x0 < 1 ) ? x0 : .5;

  for( k = 0; k < itmax; k++ ) { /* all iteration indices k */
    Int4 i;
    double g, fold = f;
    Int4 wasNewton = isNewton; /* If true, then the previous step was a */
                              /* Newton step */
    isNewton  = 0;            /* Assume that this step is not */

    /* Horner's rule for evaluating a polynomial and its derivative */
    g = 0;
    f = probs[low];
    for( i = low + d; i < 0; i += d ) {
      g = x * g + f;
      f = f * x + probs[i];
    }
    g = x * g + f;
    f = f * x + probs[0] - 1;
    for( i = d; i <= high; i += d ) {
      g = x * g + f;
      f = f * x + probs[i];
    }
    /* End Horner's rule */

    if( f > 0 ) {
      a = x; /* move the left endpoint */
    } else if( f < 0 ) {
      b = x; /* move the right endpoint */
    } else { /* f == 0 */
      break; /* x is an exact solution */
    }
    if( b - a < 2 * a * ( 1 - b ) * tolx ) {
      /* The midpoint of the interval converged */
      x = (a + b) / 2; break;
    }

    if( k >= maxNewton ||
        /* If convergence of Newton's method appears to be failing; or */
            ( wasNewton && fabs( f ) > .9 * fabs(fold) ) ||
        /* if the previous iteration was a Newton step but didn't decrease
         * f sufficiently; or */
        g >= 0
        /* if a Newton step will move us away from the desired solution */
        ) { /* then */
      /* bisect */
      x = (a + b)/2;
    } else {
      /* try a Newton step */
      double p = - f/g;
      double y = x + p;
      if( y <= a || y >= b ) { /* The proposed iterate is not in (a,b) */
        x = (a + b)/2;
      } else { /* The proposed iterate is in (a,b). Accept it. */
        isNewton = 1;
        x = y;
        if( fabs( p ) < tolx * x * (1-x) ) break; /* Converged */
      } /* else the proposed iterate is in (a,b) */
    } /* else try a Newton step. */
  } /* end for all iteration indices k */
   *itn = k;
  return -log(x)/d;
}

double UngappedBlastnKarlinBlock::KarlinLambdaNR(Blast_ScoreFreq* sfp, double initialLambdaGuess)
{
   Int4  low;        /* Lowest score (must be negative)  */
   Int4  high;       /* Highest score (must be positive) */
   Int4     itn;
   Int4  i, d;
   double*  sprob;
   double   returnValue;

   low = sfp->obs_min;
   high = sfp->obs_max;
   if (sfp->score_avg >= 0.) {   /* Expected score must be negative */
      return -1.0;
   }
   if (BlastScoreChk(low, high) != 0) return -1.;

   sprob = sfp->sprob;
   /* Find greatest common divisor of all scores */
   for (i = 1, d = -low; i <= high-low && d > 1; ++i) {
      if (sprob[i+low] != 0.0) {
         d = Gcd(d, i);
      }
   }
   returnValue =
      NlmKarlinLambdaNR( sprob, d, low, high,
                           initialLambdaGuess,
                           BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
                     20, 20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT, &itn );


   return returnValue;
}

double UngappedBlastnKarlinBlock::KarlinLtoH(Blast_ScoreFreq* sfp, double lambda)
{
   Int4  score;
   double   H, etonlam, sum, scale;

   double *probs = sfp->sprob;
   Int4 low   = sfp->obs_min,  high  = sfp->obs_max;

   if (lambda < 0.) {
      return -1.;
   }
   if (BlastScoreChk(low, high) != 0) return -1.;

   etonlam = exp( - lambda );
  sum = low * probs[low];
  for( score = low + 1; score <= high; score++ ) {
    sum = score * probs[score] + etonlam * sum;
  }

  scale = Powi( etonlam, high );
  if( scale > 0.0 ) {
    H = lambda * sum/scale;
  } else { /* Underflow of exp( -lambda * high ) */
    H = lambda * exp( lambda * high + log(sum) );
  }
   return H;
}

static double
Expm1(double x)
{
  double	absx = ABS(x);

  if (absx > .33)
    return exp(x) - 1.;

  if (absx < 1.e-16)
    return x;

  return x * (1. + x *
             (1./2. + x *
             (1./6. + x *
             (1./24. + x *
             (1./120. + x *
             (1./720. + x *
             (1./5040. + x *
             (1./40320. + x *
             (1./362880. + x *
             (1./3628800. + x *
             (1./39916800. + x *
             (1./479001600. +
              x/6227020800.))))))))))));    
}

double UngappedBlastnKarlinBlock::KarlinLHtoK(Blast_ScoreFreq* sfp, double lambda, double H)
{
    /*The next array stores the probabilities of getting each possible
      score in an alignment of fixed length; the array is shifted
      during part of the computation, so that
      entry 0 is for score 0.  */
    double         *alignmentScoreProbabilities = NULL;
    Int4            low;    /* Lowest score (must be negative) */
    Int4            high;   /* Highest score (must be positive) */
    Int4            range;  /* range of scores, computed as high - low*/
    double          K;      /* local copy of K  to return*/
    int             i;   /*loop index*/
    int             iterCounter; /*counter on iterations*/
    Int4            divisor; /*candidate divisor of all scores with
                               non-zero probabilities*/
    /*highest and lowest possible alignment scores for current length*/
    Int4            lowAlignmentScore, highAlignmentScore;
    Int4            first, last; /*loop indices for dynamic program*/
    register double innerSum;
    double          oldsum, oldsum2;  /* values of innerSum on previous
                                         iterations*/
    double          outerSum;        /* holds sum over j of (innerSum
                                        for iteration j/j)*/

    double          score_avg; /*average score*/
    /*first term to use in the closed form for the case where
      high == 1 or low == -1, but not both*/
    double          firstTermClosedForm;  /*usually store H/lambda*/
    int             iterlimit; /*upper limit on iterations*/
    double          sumlimit; /*lower limit on contributions
                                to sum over scores*/

    /*array of score probabilities reindexed so that low is at index 0*/
    double         *probArrayStartLow;

    /*pointers used in dynamic program*/
    double         *ptrP, *ptr1, *ptr2, *ptr1e;
    double          expMinusLambda; /*e^^(-Lambda) */

    if (lambda <= 0. || H <= 0.) {
        /* Theory dictates that H and lambda must be positive, so
         * return -1 to indicate an error */
        return -1.;
    }

    /*Karlin-Altschul theory works only if the expected score
      is negative*/
    if (sfp->score_avg >= 0.0) {
        return -1.;
    }

    low   = sfp->obs_min;
    high  = sfp->obs_max;
    range = high - low;

    probArrayStartLow = &sfp->sprob[low];
    /* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
       Karlin&Altschul (1990) */
    for (i = 1, divisor = -low; i <= range && divisor > 1; ++i) {
        if (probArrayStartLow[i] != 0.0)
            divisor = Gcd(divisor, i);
    }

    high   /= divisor;
    low    /= divisor;
    lambda *= divisor;

    range = high - low;

    firstTermClosedForm = H/lambda;
    expMinusLambda      = exp((double) -lambda);

    if (low == -1 && high == 1) {
        K = (sfp->sprob[low*divisor] - sfp->sprob[high*divisor]) *
            (sfp->sprob[low*divisor] - sfp->sprob[high*divisor]) / sfp->sprob[low*divisor];
        return(K);
    }

    if (low == -1 || high == 1) {
        if (high != 1) {
            score_avg = sfp->score_avg / divisor;
            firstTermClosedForm
                = (score_avg * score_avg) / firstTermClosedForm;
        }
        return firstTermClosedForm * (1.0 - expMinusLambda);
    }

    sumlimit  = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
    iterlimit = BLAST_KARLIN_K_ITER_MAX;

    alignmentScoreProbabilities =
        (double *)calloc((iterlimit*range + 1), sizeof(*alignmentScoreProbabilities));
    if (alignmentScoreProbabilities == NULL)
        return -1.;

    outerSum = 0.;
    lowAlignmentScore = highAlignmentScore = 0;
    alignmentScoreProbabilities[0] = innerSum = oldsum = oldsum2 = 1.;

    for (iterCounter = 0;
         ((iterCounter < iterlimit) && (innerSum > sumlimit));
         outerSum += innerSum /= ++iterCounter) {
        first = last = range;
        lowAlignmentScore  += low;
        highAlignmentScore += high;
        /*dynamic program to compute P(i,j)*/
        for (ptrP = alignmentScoreProbabilities +
                 (highAlignmentScore-lowAlignmentScore);
             ptrP >= alignmentScoreProbabilities;
             *ptrP-- =innerSum) {
            ptr1  = ptrP - first;
            ptr1e = ptrP - last;
            ptr2  = probArrayStartLow + first;
            for (innerSum = 0.; ptr1 >= ptr1e; ) {
                innerSum += *ptr1  *  *ptr2;
		ptr1--;
		ptr2++;
            }
            if (first)
                --first;
            if (ptrP - alignmentScoreProbabilities <= range)
                --last;
        }
        /* Horner's rule */
        innerSum = *++ptrP;
        for( i = lowAlignmentScore + 1; i < 0; i++ ) {
            innerSum = *++ptrP + innerSum * expMinusLambda;
        }
        innerSum *= expMinusLambda;

        for (; i <= highAlignmentScore; ++i)
            innerSum += *++ptrP;
        oldsum2 = oldsum;
        oldsum  = innerSum;
    }

#ifdef ADD_GEOMETRIC_TERMS_TO_K
    /*old code assumed that the later terms in sum were
      asymptotically comparable to those of a geometric
      progression, and tried to speed up convergence by
      guessing the estimated ratio between sucessive terms
      and using the explicit terms of a geometric progression
      to speed up convergence. However, the assumption does not
      always hold, and convergenece of the above code is fast
      enough in practice*/
    /* Terms of geometric progression added for correction */
    {
        double     ratio;  /* fraction used to generate the
                                   geometric progression */

        ratio = oldsum / oldsum2;
        if (ratio >= (1.0 - sumlimit*0.001)) {
            K = -1.;
            if (alignmentScoreProbabilities != NULL)
                sfree(alignmentScoreProbabilities);
            return K;
        }
        sumlimit *= 0.01;
        while (innerSum > sumlimit) {
            oldsum   *= ratio;
            outerSum += innerSum = oldsum / ++iterCounter;
        }
    }
#endif

    K = -exp((double)-2.0*outerSum) /
             (firstTermClosedForm * Expm1(-(double)lambda));

    if (alignmentScoreProbabilities != NULL)
        {free(alignmentScoreProbabilities); alignmentScoreProbabilities = NULL; }

    return K;
}

Int2 UngappedBlastnKarlinBlock::KarlinBlkUngappedCalc(Blast_KarlinBlk* kbp, Blast_ScoreFreq* sfp)
{
   if (kbp == NULL || sfp == NULL)
      return 1;

   /* Calculate the parameter Lambda */

   kbp->Lambda = KarlinLambdaNR(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT);
   if (kbp->Lambda < 0.)
      goto ErrExit;


   /* Calculate H */

   kbp->H = KarlinLtoH(sfp, kbp->Lambda);
   if (kbp->H < 0.)
      goto ErrExit;


   /* Calculate K and log(K) */

   kbp->K = KarlinLHtoK(sfp, kbp->Lambda, kbp->H);
   if (kbp->K < 0.)
      goto ErrExit;
   kbp->logK = log(kbp->K);

   /* Normal return */
   return 0;

ErrExit:
   kbp->Lambda = kbp->H = kbp->K = -1.;
   kbp->logK = HUGE_VAL;
   return 1;
}

/* gapped_karlin_block_functions */
/** Supported substitution and gap costs with corresponding quality values
 * for nucleotide sequence comparisons.
 * NB: the values 0 and 0 for the gap costs are treated as the defaults used for
 * the greedy gapped extension, i.e.
 * gap opening = 0,
 * gap extension = 1/2 match - mismatch.
 *
 * The fields are:
 *
 * 1. Gap opening cost,
 * 2. Gap extension cost,
 * 3. Lambda,
 * 4. K,
 * 5. H,
 * 6. Alpha,
 * 7. Beta,
 * 8. Theta
 */

// TODO: add gumbel parameters for nucleotide cases

/** Karlin-Altschul parameter values for substitution scores 1 and -5. */
static const array_of_8 blastn_values_1_5[] = {
    { 0, 0, 1.39, 0.747, 1.38, 1.00,  0, 100 },
    { 3, 3, 1.39, 0.747, 1.38, 1.00,  0, 100 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -4. */
static const array_of_8 blastn_values_1_4[] = {
    { 0, 0, 1.383, 0.738, 1.36, 1.02,  0, 100 },
    { 1, 2,  1.36,  0.67,  1.2,  1.1,  0,  98 },
    { 0, 2,  1.26,  0.43, 0.90,  1.4, -1,  91 },
    { 2, 1,  1.35,  0.61,  1.1,  1.2, -1,  98 },
    { 1, 1,  1.22,  0.35, 0.72,  1.7, -3,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -7.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
static const array_of_8 blastn_values_2_7[] = {
    { 0, 0,  0.69, 0.73, 1.34, 0.515,  0, 100 },
    { 2, 4,  0.68, 0.67,  1.2,  0.55,  0,  99 },
    { 0, 4,  0.63, 0.43, 0.90,   0.7, -1,  91 },
    { 4, 2, 0.675, 0.62,  1.1,   0.6, -1,  98 },
    { 2, 2,  0.61, 0.35, 0.72,   1.7, -3,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -3. */
static const array_of_8 blastn_values_1_3[] = {
    { 0, 0, 1.374, 0.711, 1.31, 1.05,  0, 100 },
    { 2, 2,  1.37,  0.70,  1.2,  1.1,  0,  99 },
    { 1, 2,  1.35,  0.64,  1.1,  1.2, -1,  98 },
    { 0, 2,  1.25,  0.42, 0.83,  1.5, -2,  91 },
    { 2, 1,  1.34,  0.60,  1.1,  1.2, -1,  97 },
    { 1, 1,  1.21,  0.34, 0.71,  1.7, -2,  88 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -5.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
static const array_of_8 blastn_values_2_5[] = {
    { 0, 0, 0.675, 0.65,  1.1,  0.6, -1, 99 },
    { 2, 4,  0.67, 0.59,  1.1,  0.6, -1, 98 },
    { 0, 4,  0.62, 0.39, 0.78,  0.8, -2, 91 },
    { 4, 2,  0.67, 0.61,  1.0, 0.65, -2, 98 },
    { 2, 2,  0.56, 0.32, 0.59, 0.95, -4, 82 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -2. */
static const array_of_8 blastn_values_1_2[] = {
    { 0, 0, 1.28, 0.46, 0.85, 1.5, -2, 96 },
    { 2, 2, 1.33, 0.62,  1.1, 1.2,  0, 99 },
    { 1, 2, 1.30, 0.52, 0.93, 1.4, -2, 97 },
    { 0, 2, 1.19, 0.34, 0.66, 1.8, -3, 89 },
    { 3, 1, 1.32, 0.57,  1.0, 1.3, -1, 99 },
    { 2, 1, 1.29, 0.49, 0.92, 1.4, -1, 96 },
    { 1, 1, 1.14, 0.26, 0.52, 2.2, -5, 85 }
};

/** Karlin-Altschul parameter values for substitution scores 2 and -3.
 * These parameters can only be applied to even scores. Any odd score must be
 * rounded down to the nearest even number before calculating the e-value.
 */
static const array_of_8 blastn_values_2_3[] = {
    { 0, 0,  0.55, 0.21, 0.46,  1.2, -5, 87 },
    { 4, 4,  0.63, 0.42, 0.84, 0.75, -2, 99 },
    { 2, 4, 0.615, 0.37, 0.72, 0.85, -3, 97 },
    { 0, 4,  0.55, 0.21, 0.46,  1.2, -5, 87 },
    { 3, 3, 0.615, 0.37, 0.68,  0.9, -3, 97 },
    { 6, 2,  0.63, 0.42, 0.84, 0.75, -2, 99 },
    { 5, 2, 0.625, 0.41, 0.78,  0.8, -2, 99 },
    { 4, 2,  0.61, 0.35, 0.68,  0.9, -3, 96 },
    { 2, 2, 0.515, 0.14, 0.33, 1.55, -9, 81 }
};

/** Karlin-Altschul parameter values for substitution scores 3 and -4. */
static const array_of_8 blastn_values_3_4[] = {
    { 6, 3, 0.389, 0.25, 0.56, 0.7, -5, 95},
    { 5, 3, 0.375, 0.21, 0.47, 0.8, -6, 92},
    { 4, 3, 0.351, 0.14, 0.35, 1.0, -9, 86},
    { 6, 2, 0.362, 0.16, 0.45, 0.8, -4, 88},
    { 5, 2, 0.330, 0.092, 0.28, 1.2, -13, 81},
    { 4, 2, 0.281, 0.046, 0.16, 1.8, -23, 69}
};

/** Karlin-Altschul parameter values for substitution scores 4 and -5. */
static const array_of_8 blastn_values_4_5[] = {
    { 0, 0, 0.22, 0.061, 0.22, 1.0, -15, 74 },
    { 6, 5, 0.28,  0.21, 0.47, 0.6 , -7, 93 },
    { 5, 5, 0.27,  0.17, 0.39, 0.7,  -9, 90 },
    { 4, 5, 0.25,  0.10, 0.31, 0.8, -10, 83 },
    { 3, 5, 0.23, 0.065, 0.25, 0.9, -11, 76 }
};

/** Karlin-Altschul parameter values for substitution scores 1 and -1. */
static const array_of_8 blastn_values_1_1[] = {
    { 3,  2, 1.09,  0.31, 0.55, 2.0,  -2, 99 },
    { 2,  2, 1.07,  0.27, 0.49, 2.2,  -3, 97 },
    { 1,  2, 1.02,  0.21, 0.36, 2.8,  -6, 92 },
    { 0,  2, 0.80, 0.064, 0.17, 4.8, -16, 72 },
    { 4,  1, 1.08,  0.28, 0.54, 2.0,  -2, 98 },
    { 3,  1, 1.06,  0.25, 0.46, 2.3,  -4, 96 },
    { 2,  1, 0.99,  0.17, 0.30, 3.3, -10, 90 }
};

/** Karlin-Altschul parameter values for substitution scores 3 and -2. */
static const array_of_8 blastn_values_3_2[] = {
    {  5,  5, 0.208, 0.030, 0.072, 2.9, -47, 77}
};

/** Karlin-Altschul parameter values for substitution scores 5 and -4. */
static const array_of_8 blastn_values_5_4[] = {
    { 10, 6, 0.163, 0.068, 0.16, 1.0, -19, 85 },
    {  8, 6, 0.146, 0.039, 0.11, 1.3, -29, 76 }
};

Int2 GappedBlastnKarlinBlock::SplitArrayOf8(const array_of_8* input, const array_of_8** normal,
                                            const array_of_8** non_affine, Boolean* split)
{
    if (input == NULL || normal == NULL || non_affine == NULL)
        return -1;

    *normal = NULL;
    *non_affine = NULL;

    // gap_open == 0 && gap_extend == 0
    if (input[0][0] == 0 && input[0][1] == 0)
    {
        *normal = input + 1;
        *non_affine = input;
        *split = TRUE;
    }
    else
    {
        *normal = input;
        *split = FALSE;
    }
    return 0;
}

Int2 GappedBlastnKarlinBlock::AdjustGapParametersByGcd(array_of_8* normal, array_of_8* linear,
            int size, Int4* gap_existence_max, Int4* gap_extend_max, int divisor)
{
    if (divisor == 1) return 0;

    if (size <= 0) return 1;

    (*gap_existence_max) *= divisor;
    (*gap_extend_max) *= divisor;

    if (normal)
    {
        int i;
        for (i = 0; i < size; ++i)
        {
            normal[i][0] *= divisor;
            normal[i][1] *= divisor;
            normal[i][2] /= divisor;
            normal[i][5] /= divisor;
        }
    }
    if (linear)
    {
        linear[0][0] *= divisor;
        linear[0][1] *= divisor;
        linear[0][2] /= divisor;
        linear[0][5] /= divisor;
    }

    return 0;
}

Int2 GappedBlastnKarlinBlock::GetNuclValuesArray(Int4 reward, Int4 penalty,
                                                 Int4* array_size, array_of_8** normal,
                                                 array_of_8** non_affine, Int4* gap_open_max,
                                                 Int4* gap_extend_max, Boolean* round_down)
{
    Int2 status = 0;
    const array_of_8* kValues = NULL;
    const array_of_8* kValues_non_affine = NULL;
    Boolean split = FALSE;
    int divisor = Gcd(reward, penalty);

    *round_down = FALSE;

    *array_size = 0;
    *normal = NULL;
    *non_affine = NULL;

    if (divisor != 0)
    {
        reward /= divisor;
        penalty /= divisor;
    }

    if (reward == 1 && penalty == -5) {
        if ((status=SplitArrayOf8(blastn_values_1_5, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_1_5)/sizeof(array_of_8);
        *gap_open_max = 3;
        *gap_extend_max = 3;
    } else if (reward == 1 && penalty == -4) {
        if ((status=SplitArrayOf8(blastn_values_1_4, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_1_4)/sizeof(array_of_8);
        *gap_open_max = 2;
        *gap_extend_max = 2;
    } else if (reward == 2 && penalty == -7) {
        if ((status=SplitArrayOf8(blastn_values_2_7, &kValues, &kValues_non_affine, &split)))
           return status;

        *round_down = TRUE;
        *array_size = sizeof(blastn_values_2_7)/sizeof(array_of_8);
        *gap_open_max = 4;
        *gap_extend_max = 4;
    } else if (reward == 1 && penalty == -3) {
        if ((status=SplitArrayOf8(blastn_values_1_3, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_1_3)/sizeof(array_of_8);
        *gap_open_max = 2;
        *gap_extend_max = 2;
    } else if (reward == 2 && penalty == -5) {
        if ((status=SplitArrayOf8(blastn_values_2_5, &kValues, &kValues_non_affine, &split)))
           return status;

        *round_down = TRUE;
        *array_size = sizeof(blastn_values_2_5)/sizeof(array_of_8);
        *gap_open_max = 4;
        *gap_extend_max = 4;
    } else if (reward == 1 && penalty == -2) {
        if ((status=SplitArrayOf8(blastn_values_1_2, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_1_2)/sizeof(array_of_8);
        *gap_open_max = 2;
        *gap_extend_max = 2;
    } else if (reward == 2 && penalty == -3) {
        if ((status=SplitArrayOf8(blastn_values_2_3, &kValues, &kValues_non_affine, &split)))
           return status;

        *round_down = TRUE;
        *array_size = sizeof(blastn_values_2_3)/sizeof(array_of_8);
        *gap_open_max = 6;
        *gap_extend_max = 4;
    } else if (reward == 3 && penalty == -4) {
        if ((status=SplitArrayOf8(blastn_values_3_4, &kValues, &kValues_non_affine, &split)))
           return status;

        *round_down = TRUE;
        *array_size = sizeof(blastn_values_3_4)/sizeof(array_of_8);
        *gap_open_max = 6;
        *gap_extend_max = 3;
    } else if (reward == 1 && penalty == -1) {
        if ((status=SplitArrayOf8(blastn_values_1_1, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_1_1)/sizeof(array_of_8);
        *gap_open_max = 4;
        *gap_extend_max = 2;
    } else if (reward == 3 && penalty == -2) {
        if ((status=SplitArrayOf8(blastn_values_3_2, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_3_2)/sizeof(array_of_8);
        *gap_open_max = 5;
        *gap_extend_max = 5;
    } else if (reward == 4 && penalty == -5) {
        if ((status=SplitArrayOf8(blastn_values_4_5, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_4_5)/sizeof(array_of_8);
        *gap_open_max = 12;
        *gap_extend_max = 8;
    } else if (reward == 5 && penalty == -4) {
        if ((status=SplitArrayOf8(blastn_values_5_4, &kValues, &kValues_non_affine, &split)))
           return status;

        *array_size = sizeof(blastn_values_5_4)/sizeof(array_of_8);
        *gap_open_max = 25;
        *gap_extend_max = 10;
    } else  { /* Unsupported reward-penalty */
        status = -1;
        fprintf(stderr, "[%s] Substitution scores %d and %d are not supported.\n", __func__, reward, penalty);
        exit(1);
    }

    if (split)
        (*array_size)--;

    if (status == 0)
    {
        if (*array_size > 0)
        {
            void* copy = malloc((*array_size) * sizeof(array_of_8));
            ASSERT(copy != NULL);
            memcpy(copy, kValues, (*array_size) * sizeof(array_of_8));
            *normal = (array_of_8*)copy;
        }
        if (kValues_non_affine)
        {
            void* copy = malloc(sizeof(array_of_8));
            ASSERT(copy != NULL);
            memcpy(copy, kValues_non_affine, sizeof(array_of_8));
            *non_affine = (array_of_8*)copy;
        }

        status = AdjustGapParametersByGcd(*normal, *non_affine, *array_size, gap_open_max, gap_extend_max, divisor);
    }

    return status;
}

Int2 GappedBlastnKarlinBlock::KarlinBlkNuclGappedCalc(Blast_KarlinBlk* kbp,
                                                      Int4 gap_open, Int4 gap_extend,
                                                      Int4 reward, Int4 penalty,
                                                      Blast_KarlinBlk* kbp_ungap,
                                                      Boolean* round_down)
{
    const int kGapOpenIndex = 0;
    const int kGapExtIndex = 1;
    const int kLambdaIndex = 2;
    const int kKIndex = 3;
    const int kHIndex = 4;
    int num_combinations = 0;
    int gap_open_max, gap_extend_max;
    array_of_8* normal=NULL;
    array_of_8* linear=NULL;
    Int2 status  =  GetNuclValuesArray(reward,
                                       penalty,
                                       &num_combinations,
                                       &normal,
                                       &linear,
                                       &gap_open_max,
                                       &gap_extend_max,
                                       round_down);

    if (status)
    {
       free(normal); normal = NULL;
       free(linear); linear = NULL;
       return status;
    }

    ASSERT(kbp && kbp_ungap);


    /* Try to find the table entry corresponding to input gap costs. */
    if (gap_open == 0 && gap_extend == 0 && linear)
    {
        kbp->Lambda = linear[0][kLambdaIndex];
        kbp->K = linear[0][kKIndex];
        kbp->logK = log(kbp->K);
        kbp->H = linear[0][kHIndex];
    }
    else
    {
        int index=0;
        for (index = 0; index < num_combinations; ++index) {
            if (normal[index][kGapOpenIndex] == gap_open &&
                normal[index][kGapExtIndex] == gap_extend) {
                kbp->Lambda = normal[index][kLambdaIndex];
                kbp->K = normal[index][kKIndex];
                kbp->logK = log(kbp->K);
                kbp->H = normal[index][kHIndex];
                break;
            }
        }

        /* If gap costs are not found in the table, check if they belong to the
        infinite domain, where ungapped values of the parameters can be used. */
        if (index == num_combinations) {
        /* If gap costs are larger than maximal provided in tables, copy
           the values from the ungapped Karlin block. */
            if (gap_open >= gap_open_max && gap_extend >= gap_extend_max) {
                Blast_KarlinBlkCopy(kbp, kbp_ungap);
            } else {
                int i = 0;
                /* Unsupported gap costs combination. */
                fprintf(stderr, "Gap existence and extension values %ld and %ld "
                        "are not supported for substitution scores %ld and %ld\n",
                        (long) gap_open, (long) gap_extend, (long) reward, (long) penalty);
                for (i = 0; i < num_combinations; ++i)
                {
                     fprintf(stderr, "%ld and %ld are supported existence and extension values\n",
                        (long) normal[i][kGapOpenIndex],  (long) normal[i][kGapExtIndex]);
                }
                fprintf(stderr, "%ld and %ld are supported existence and extension values\n",
                     (long) gap_open_max, (long) gap_extend_max);
                fprintf(stderr, "Any values more stringent than %ld and %ld are supported\n",
                     (long) gap_open_max, (long) gap_extend_max);
                free(normal); normal = NULL;
                free(linear); linear = NULL;
                return 1;
            }
        }
    }

    free(normal); normal = NULL;
    free(linear); linear = NULL;
    return 0;
}

/* BlastScoreBlk functions */

BlastScoreBlk::BlastScoreBlk(Uint1 alphabet, QueryInfo* queries)
{
    alphabet_code = alphabet;
    alphabet_size = BLASTNA_SIZE;
    alphabet_start = 0;

    ASSERT(queries != NULL);
    query_blk = queries;
    round_down = FALSE;

    matrix = SBlastScoreMatrixNew(alphabet_size, alphabet_size);
    ASSERT(matrix != NULL);

    scale_factor = 1.0;
    number_of_contexts = queries->Size() * 2;

    sfp = (Blast_ScoreFreq**)calloc(number_of_contexts, sizeof(Blast_ScoreFreq*));
    kbp = (Blast_KarlinBlk**)calloc(number_of_contexts, sizeof(Blast_KarlinBlk*));
    kbp_gap = (Blast_KarlinBlk**)calloc(number_of_contexts, sizeof(Blast_KarlinBlk*));

    eff_searchsp = (Int8*)malloc(number_of_contexts * sizeof(Int8));
    length_adjustment = (Int4*)malloc(number_of_contexts * sizeof(Int4));

    ambiguous_res = NULL;
    ambig_size = 0;
    ambig_occupy = 0;
}

BlastScoreBlk::~BlastScoreBlk()
{
    ScoreBlkFree();
}

Int2 BlastScoreBlk::ScoreBlkMatrixFill()
{
    Int2 status = 0;
    status = BlastScoreBlkNuclMatrixCreate(this);
    if (status != 0) return status;

    status = BlastScoreBlkMaxScoreSet(this);

    return status;
}

Int2 BlastScoreBlk::ScoreBlkMatrixInit(const BlastScoringOptions* scoring_options)
{
    Int2 status = 0;

    if (!scoring_options) return -1;

    ScoreBlkSetAmbigRes('N');
    ScoreBlkSetAmbigRes('-');

    penalty = scoring_options->penalty;
    reward  = scoring_options->reward;

    status = ScoreBlkMatrixFill();
    return status;
}


Int2 BlastScoreBlk::ScoreBlkInit(const BlastScoringOptions* scoring_options,
				 double sf)
{
    Int2 status = 0;

    scale_factor = sf;

    status = ScoreBlkMatrixInit(scoring_options);
    if (status)
    {
	fprintf(stderr, "[%s] Cannot initialize matrix.\n", __func__);
	return status;
    }

    status = ScoreBlkKbpUngappedCalc();

    if (scoring_options->gapped_calculation)
    {
	status = ScoreBlkKbpGappedCalc(scoring_options);
    }

    return status;
}

void BlastScoreBlk::ScoreBlkFree()
{
    Int4 index;
    for (index = 0; index < number_of_contexts; ++index)
    {
        if (sfp)
            sfp[index] = Blast_ScoreFreqFree(sfp[index]);
        if (kbp)
            kbp[index] = Blast_KarlinBlkFree(kbp[index]);
        if (kbp_gap)
            kbp_gap[index] = Blast_KarlinBlkFree(kbp_gap[index]);
    }
    if (kbp_ideal)
	kbp_ideal = Blast_KarlinBlkFree(kbp_ideal);

    sfree(sfp);
    sfree(kbp);
    sfree(kbp_gap);

    matrix = SBlastScoreMatrixFree(matrix);
    sfree(ambiguous_res);

    sfree(eff_searchsp);
    sfree(length_adjustment);
}

Int2 BlastScoreBlk::ScoreBlkKbpIdealCalc()
{
    Blast_ResFreq* stdrfp = NULL;
    Blast_ScoreFreq* sfp = NULL;
    Int2 status = 0;

    stdrfp = Blast_ResFreqNew(this);
    Blast_ResFreqStdComp(this, stdrfp);
    sfp = Blast_ScoreFreqNew(loscore, hiscore);
    BlastScoreFreqCalc(this, sfp, stdrfp, stdrfp);
    kbp_ideal = Blast_KarlinBlkNew();

    UngappedBlastnKarlinBlock ungapped_karlinblk;
    ungapped_karlinblk.KarlinBlkUngappedCalc(kbp_ideal, sfp);

    stdrfp = Blast_ResFreqFree(stdrfp);
    sfp = Blast_ScoreFreqFree(sfp);
    return status;
}
Int2 BlastScoreBlk::ScoreBlkKbpGappedCalc(const BlastScoringOptions* scoring_options)
{
    Int4 index = 0;
    Int2 status = 0;

    if (scoring_options == NULL) return 1;

    GappedBlastnKarlinBlock gapped_karlinblk;

    Int4 num_queries = query_blk->Size();
    Int4 i;
    for (index = 0; index < num_queries; ++index)
    {
	for (i = 0; i < 2; ++i)
	{
	    Int4 j = 2 * index + i;
		if (query_blk->contexts[j].is_valid == FALSE) continue;

	    if (kbp[j] == NULL) continue;

	    kbp_gap[j] = Blast_KarlinBlkNew();

	    status = gapped_karlinblk.KarlinBlkNuclGappedCalc(kbp_gap[j],
							      scoring_options->gap_open,
							      scoring_options->gap_extend,
						              scoring_options->reward,
							      scoring_options->penalty,
							      kbp[j], &round_down);
	    if (status) exit(1); //return status;
	}
    }

    return 0;
}

Int2 BlastScoreBlk::ScoreBlkKbpUngappedCalc()
{
    Int2 status = 0;
    Blast_ResFreq* rfp, *stdrfp;

    // Ideal Karlin block is filled unconditionally
    status = ScoreBlkKbpIdealCalc();
    if (status) return status;

    stdrfp = Blast_ResFreqNew(this);
    Blast_ResFreqStdComp(this, stdrfp);
    rfp = Blast_ResFreqNew(this);

    Int4 num_queries = query_blk->Size();
    Int4 i, j;
    Int4 context_offset;
    Int4 context;
    Boolean valid_context = FALSE;

    UngappedBlastnKarlinBlock ungapped_karlinblk;

    for (i = 0; i < num_queries; ++i)
    {
	for (j = 0; j < 2; ++j)
	{
	    Int4 query_length;
	    const Uint1* buffer;  // holds sequence
	    Blast_KarlinBlk* kbp_tmp;
	    Int2 loop_status;  // status flag for functions in this loop

	    query_length = query_blk->GetSeqLength(i * 2);
	    context_offset = 0;

            buffer = query_blk->GetSequence(2 * i + j);

	    Blast_ResFreqString(this, rfp, (char*)buffer, query_length);

	    context = 2 * i + j;
	    sfp[context] = Blast_ScoreFreqNew(loscore, hiscore);
	    ASSERT(sfp[context] != NULL);
	    BlastScoreFreqCalc(this, sfp[context], rfp, stdrfp);
	    kbp[context] = kbp_tmp = Blast_KarlinBlkNew();
	    loop_status = ungapped_karlinblk.KarlinBlkUngappedCalc(kbp_tmp, sfp[context]);

	    if (loop_status)
	    {
			query_blk->contexts[context].is_valid = FALSE;
			fprintf(stderr, "[%s] Warning: Could not calculate ungapped Karlin-Altschul parameters due "
							"to an invalid query sequence. Please verify the "
							"query sequence(s) and/or filtering options.\n", __func__);

			sfp[context] = Blast_ScoreFreqFree(sfp[context]);
			kbp[context] = Blast_KarlinBlkFree(kbp[context]);

			continue;
	    }
	    valid_context = TRUE;
	}
    }

    rfp = Blast_ResFreqFree(rfp);
    stdrfp = Blast_ResFreqFree(stdrfp);

    if (valid_context == FALSE)
    {
	status = 1;
    }
    return status;
}

Int2 BlastScoreBlk::ScoreBlkSetAmbigRes(char ambig_res)
{
	Int2 index;
	Uint1* ambig_buffer;

	if (ambig_occupy >= ambig_size)
	{
		ambig_size += 5;
		ambig_buffer = (Uint1*)calloc(ambig_size, sizeof(Uint1));
		memcpy(ambig_buffer, ambiguous_res, ambig_occupy * sizeof(Uint1));
		sfree(ambiguous_res);
		ambiguous_res = ambig_buffer;
	}

	if (alphabet_code == BLASTNA_SEQ_CODE)
		ambiguous_res[ambig_occupy] =
				IUPACNA_TO_BLASTNA[toupper((unsigned char) ambig_res)];
	else
		ambiguous_res[ambig_occupy] =
				IUPACNA_TO_NCBI4NA[toupper((unsigned char) ambig_res)];

	++ambig_occupy;

	return 0;
}

static Boolean
BlastEffectiveLengthOptions_IsSearchSpaceSet(
		const BlastEffectiveLengthsOptions* options)
{
	int i;
	if (!options || options->searchsp_eff == NULL)
		return FALSE;

	for (i = 0; i < options->num_searchspaces; ++i)
	{
		if (options->searchsp_eff[i] != 0)
			return TRUE;
	}

	return FALSE;
}

Int8 BlastScoreBlk::s_GetEffectiveSearchSpaceForContext(
		const BlastEffectiveLengthsOptions* eff_len_options,
		int context_index)
{
	Int8 retval = 0;

	if (eff_len_options->num_searchspaces == 0)
		retval = 0;
	else if (eff_len_options->num_searchspaces == 1)
	{
		//if (context_index != 0)
		//fprintf(stderr, "[%s] Warning: One search space is being used for multiple sequences.\n",
		//		__func__);
		//retval = eff_len_options->searchsp_eff[0];
	}
	else
	{
		//ASSERT(context_index < eff_len_options->num_searchspaces);
            if (context_index >= eff_len_options->num_searchspaces)
            {
                fprintf(stderr, "[%s] context_index = %d\n", __func__, context_index);
                fprintf(stderr, "\tnum_searchspace = %d\n", eff_len_options->num_searchspaces);
                abort();
            }
		retval = eff_len_options->searchsp_eff[context_index];
	}

	return retval;
}

/**
 * Computes the adjustment to the lengths of the query and database sequences
 * that is used to compensate for edge effects when computing evalues.
 *
 * The length adjustment is an integer-valued approximation to the fixed
 * point of the function
 *
 *    f(ell) = beta +
 *               (alpha/lambda) * (log K + log((m - ell)*(n - N ell)))
 *
 * where m is the query length n is the length of the database and N is the
 * number of sequences in the database. The values beta, alpha, lambda and
 * K are statistical, Karlin-Altschul parameters.
 *
 * The value of the length adjustment computed by this routine, A,
 * will always be an integer smaller than the fixed point of
 * f(ell). Usually, it will be the largest such integer.  However, the
 * computed length adjustment, A, will also be so small that
 *
 *    K * (m - A) * (n - N * A) > MAX(m,n).
 *
 * Moreover, an iterative method is used to compute A, and under
 * unusual circumstances the iterative method may not converge.
 *
 * @param K      the statistical parameter K
 * @param logK   the natural logarithm of K
 * @param alpha_d_lambda    the ratio of the statistical parameters
 *                          alpha and lambda (for ungapped alignments, the
 *                          value 1/H should be used)
 * @param beta              the statistical parameter beta (for ungapped
 *                          alignments, beta == 0)
 * @param query_length      the length of the query sequence
 * @param db_length         the length of the database
 * @param db_num_seqs       the number of sequences in the database
 * @param length_adjustment the computed value of the length adjustment [out]
 *
 * @return   0 if length_adjustment is known to be the largest integer less
 *           than the fixed point of f(ell); 1 otherwise.
 */
Int4
BLAST_ComputeLengthAdjustment(double K,
                             double logK,
                             double alpha_d_lambda,
                             double beta,
                             Int4 query_length,
                             Int8 db_length,
                             Int4 db_num_seqs,
                             Int4 * length_adjustment)
{
    Int4 i;                     /* iteration index */
    const Int4 kMaxIterations = 20;     /* maximum allowed iterations */
    double m = (double) query_length;
    double n = (double) db_length;
    double N = (double) db_num_seqs;

    double ell;                 /* A float value of the length adjustment */
    double ss;                  /* effective size of the search space */
    double ell_min = 0, ell_max;   /* At each iteration i,
                                         * ell_min <= ell <= ell_max. */
    Boolean converged    = FALSE;       /* True if the iteration converged */
    double ell_next = 0;        /* Value the variable ell takes at iteration
                                 * i + 1 */
    /* Choose ell_max to be the largest nonnegative value that satisfies
     *
     *    K * (m - ell) * (n - N * ell) > MAX(m,n)
     *
     * Use quadratic formula: 2 c /( - b + sqrt( b*b - 4 * a * c )) */
    { /* scope of a, mb, and c, the coefficients in the quadratic formula
       * (the variable mb is -b) */
        double a  = N;
        double mb = m * N + n;
        double c  = n * m - MAX(m, n) / K;

        if(c < 0) {
            *length_adjustment = 0;
            return 1;
        } else {
            ell_max = 2 * c / (mb + sqrt(mb * mb - 4 * a * c));
        }
    } /* end scope of a, mb and c */

    for(i = 1; i <= kMaxIterations; i++) {      /* for all iteration indices */
        double ell_bar;         /* proposed next value of ell */
        ell      = ell_next;
        ss       = (m - ell) * (n - N * ell);
        ell_bar  = alpha_d_lambda * (logK + log(ss)) + beta;
        if(ell_bar >= ell) { /* ell is no bigger than the true fixed point */
            ell_min = ell;
            if(ell_bar - ell_min <= 1.0) {
                converged = TRUE;
                break;
            }
            if(ell_min == ell_max) { /* There are no more points to check */
                break;
            }
        } else { /* else ell is greater than the true fixed point */
            ell_max = ell;
        }
        if(ell_min <= ell_bar && ell_bar <= ell_max) {
          /* ell_bar is in range. Accept it */
            ell_next = ell_bar;
        } else { /* else ell_bar is not in range. Reject it */
            ell_next = (i == 1) ? ell_max : (ell_min + ell_max) / 2;
        }
    } /* end for all iteration indices */
    if(converged) { /* the iteration converged */
        /* If ell_fixed is the (unknown) true fixed point, then we
         * wish to set (*length_adjustment) to floor(ell_fixed).  We
         * assume that floor(ell_min) = floor(ell_fixed) */
        *length_adjustment = (Int4) ell_min;
        /* But verify that ceil(ell_min) != floor(ell_fixed) */
        ell = ceil(ell_min);
        if( ell <= ell_max ) {
          ss = (m - ell) * (n - N * ell);
          if(alpha_d_lambda * (logK + log(ss)) + beta >= ell) {
            /* ceil(ell_min) == floor(ell_fixed) */
            *length_adjustment = (Int4) ell;
          }
        }
    } else { /* else the iteration did not converge. */
        /* Use the best value seen so far */
        *length_adjustment = (Int4) ell_min;
    }

    return converged ? 0 : 1;
}

Int2 BlastScoreBlk::BLAST_CalcEffLength(
		const BlastScoringOptions* scoring_options,
		const BlastEffectiveLengthsParameters* eff_len_params,
		QueryInfo* query_info)
{
	double alpha = 0.0, beta = 0.0;
	Int4 index;
	Int4 db_num_seqs;
	Int8 db_length;
	Blast_KarlinBlk** kbp_ptr;
	const BlastEffectiveLengthsOptions* eff_len_options =
			eff_len_params->options;

	if (!query_info)
		return -1;

	if (eff_len_options->db_length > 0)
		db_length = eff_len_options->db_length;
	else
		db_length = eff_len_params->real_db_length;

	if (db_length == 0
			&& !BlastEffectiveLengthOptions_IsSearchSpaceSet(eff_len_options))
		return 0;

	if (eff_len_options->dbseq_num > 0)
		db_num_seqs = eff_len_options->dbseq_num;
	else
		db_num_seqs = eff_len_params->real_num_seqs;

	kbp_ptr = (scoring_options->gapped_calculation ? kbp_gap : kbp);
        
        Int4 num_queries = query_info->Size();

	for (index = 0; index < num_queries * 2; ++index)
	{
		Blast_KarlinBlk* kbp_tmp;
		Int4 len_adjust = 0;
		Int4 query_length = query_info->GetSeqLength(index);

		Int8 effective_search_space = s_GetEffectiveSearchSpaceForContext(
				eff_len_options, index);

		kbp_tmp = kbp_ptr[index];

		if (query_blk->contexts[index].is_valid && kbp_tmp != NULL && query_length > 0)
		{
			Blast_GetNuclAlphaBeta(scoring_options->reward,
					scoring_options->penalty, scoring_options->gap_open,
					scoring_options->gap_extend, kbp[index],
					scoring_options->gapped_calculation, &alpha, &beta);

			BLAST_ComputeLengthAdjustment(kbp_tmp->K, kbp_tmp->logK,
					alpha / kbp_tmp->Lambda, beta, query_length, db_length,
					db_num_seqs, &len_adjust);

			if (effective_search_space == 0)
			{
				Int8 effective_db_length = db_length
						- ((Int8) db_num_seqs * len_adjust);
				if (effective_db_length <= 0)
					effective_db_length = 1;

				effective_search_space = effective_db_length
						* (query_length - len_adjust);
			}
		}
		eff_searchsp[index] = effective_search_space;
		length_adjustment[index] = len_adjust;
	}

	return 0;
}
