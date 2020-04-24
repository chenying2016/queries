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

#ifndef BLAST_STAT_FUNCTIONS_H
#define	BLAST_STAT_FUNCTIONS_H

#include "stat.h"
#include "sequence.h"
#include "options.h"

/** Allocates a new SBlastScoreMatrix structure of the specified dimensions.
 * @param ncols number of columns [in]
 * @param nrows number of rows [in]
 * @return NULL in case of memory allocation failure, else new
 * SBlastScoreMatrix structure
 */
SBlastScoreMatrix*
SBlastScoreMatrixNew(size_t ncols, size_t nrows);

SBlastScoreMatrix*
SBlastScoreMatrixFree(SBlastScoreMatrix* matrix);

/** This function fills in the BlastScoreBlk structure.
 * Tasks are:
 * -read in the matrix
 * -set maxscore
 * @param sbp Scoring block [in] [out]
 * @param get_path pointer to function that will return path to matrix.  Only called
 *  if built-in matrix not found [in]
*/
Int2 Blast_ScoreBlkMatrixFill (BlastScoreBlk* sbp);

Blast_KarlinBlk* Blast_KarlinBlkNew (void);

Int2 Blast_KarlinBlkCopy(Blast_KarlinBlk* kbp_to, Blast_KarlinBlk* kbp_from);

Blast_KarlinBlk* Blast_KarlinBlkFree(Blast_KarlinBlk* kbp);

/** Returns the first valid Karlin-Altchul block from the list of blocks.
 * @sa s_BlastKarlinBlkIsValid
 * @param kbp_in array of Karlin blocks to be searched [in]
 * @param query_info information on number of queries (specifies number of
 * elements in above array) [in]
 * @param kbp_ret the object to be pointed at [out]
 * @return zero on success, 1 if no valid block found
 */
Int2
s_BlastFindValidKarlinBlk(Blast_KarlinBlk** kbp_in,
		const QueryInfo* query_info,
		Blast_KarlinBlk** kbp_ret);

/** Returns true if the Karlin-Altschul block doesn't have its lambda, K, and H
 * fields set to negative values. -1 is the sentinel used to mark them as
 * invalid. This can happen if a query sequence is completely masked for
 * example.
 * @param kbp Karlin-Altschul block to examine [in]
 * @return TRUE if its valid, else FALSE
 */
Boolean s_BlastKarlinBlkIsValid(const Blast_KarlinBlk* kbp);

/** Returns the smallest lambda value from a collection
 *  of Karlin-Altchul blocks
 * @param kbp_in array of Karlin blocks to be searched [in]
 * @param query_info information on number of queries (specifies number of
 * elements in above array) [in]
 * @param kbp_out Karlin blocks with smallest lambda [out]
 * @return The smallest lambda value
 */
double
s_BlastFindSmallestLambda(Blast_KarlinBlk** kbp_in,
                          const QueryInfo* query_info,
                          Blast_KarlinBlk** kbp_out);

/** Compute a divisor used to weight the evalue of a collection of
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
 * interval (0,1). Typical values of decayrate are .1 and .5.
 * @param decayrate adjusts for (multiple) tests of number of HSP sum groups [in]
 * @param nsegs the number of HSPs in the sum group [in]
 * @return divisor used to compensate for multiple tests
 */
double BLAST_GapDecayDivisor(double decayrate, unsigned nsegs );

/** Calculates the Expect value based upon the search space and some Karlin-Altschul
 * parameters.  It is "simple" as it does not use sum-statistics.
 * @param S the score of the alignment. [in]
 * @param kbp the Karlin-Altschul parameters. [in]
 * @param searchsp total search space to be used [in]
 * @return the expect value
 */
double BLAST_KarlinStoE_simple (Int4 S, Blast_KarlinBlk* kbp, Int8  searchsp);

/** Calculate the cutoff score from the expected number of HSPs or vice versa.
 * @param S The calculated score [in] [out]
 * @param E The calculated e-value [in] [out]
 * @param kbp The Karlin-Altschul statistical parameters [in]
 * @param searchsp The effective search space [in]
 * @param dodecay Use gap decay feature? [in]
 * @param gap_decay_rate Gap decay rate to use, if dodecay is set [in]
 */
Int2 BLAST_Cutoffs (Int4 *S, double* E, Blast_KarlinBlk* kbp,
                    Int8 searchsp, Boolean dodecay, double gap_decay_rate);

/** Extract the recommended gap existence and extension values.
 * Only to be used with blastn searches.
 * @param reward match score [in]
 * @param penalty mismatch score [in]
 * @param gap_existence returns recommended existence cost [in|out]
 * @param gap_extension returns recommended extension cost [in|out]
 * @return zero on success
 */
Int2 BLAST_GetNucleotideGapExistenceExtendParams(Int4 reward,
                                       Int4 penalty,
                                       Int4* gap_existence,
                                       Int4* gap_extension);

/** Check the validity of the reward and penalty scores.
 * Only to be used with blastn searches.
 * @param reward match score [in]
 * @param penalty mismatch score [in]
 * @return TRUE on success
 */
Boolean BLAST_CheckRewardPenaltyScores(Int4 reward, Int4 penalty);

/** Extract the alpha and beta settings for this matrixName, and these
 *  gap open and gap extension costs
 * @param matrixName name of the matrix used [in]
 * @param alpha Karlin-Altschul parameter to be set [out]
 * @param beta Karlin-Altschul parameter to be set [out]
 * @param gapped TRUE if a gapped search [in]
 * @param gap_open existence cost of a gap [in]
 * @param gap_extend extension cost of a gap [in]
 * @param kbp_ungapped Karlin block with ungapped values of the parameters [in]
*/
void BLAST_GetAlphaBeta (const char* matrixName, double *alpha,
                         double *beta, Boolean gapped, Int4 gap_open,
                         Int4 gap_extend, const Blast_KarlinBlk* kbp_ungapped);

/** Extract the alpha and beta settings for these substitution and gap scores.
 * If substitution or gap costs are not found in the tables, assume an ungapped
 * search. Then alpha is computed using the formula Alpha = Lambda/H, and beta
 * is equal to 0 except for some special cases.
 * @param reward Match reward score [in]
 * @param penalty Mismatch penalty score [in]
 * @param gap_open Gap opening (existence) cost [in]
 * @param gap_extend Gap extension cost [in]
 * @param kbp Karlin block containing already computed Lambda, K and H
 *            parameters.
 * @param gapped_calculation Is this a gapped search? [in]
 * @param alpha Alpha parameter for this scoring system [out]
 * @param beta Beta parameter for this scoring system [out]
 */
Int2 Blast_GetNuclAlphaBeta(Int4 reward, Int4 penalty, Int4 gap_open,
                            Int4 gap_extend, Blast_KarlinBlk* kbp,
                            Boolean gapped_calculation,
                            double *alpha, double *beta);

Blast_ResFreq* Blast_ResFreqNew(const BlastScoreBlk* sbp);

Blast_ResFreq* Blast_ResFreqFree(Blast_ResFreq* rfp);

/** Calculates residues frequencies given a standard distribution.
 * @param sbp the BlastScoreBlk provides information on alphabet.
 * @param rfp the prob element on this Blast_ResFreq is used.
 * @return zero on success
*/
Int2 Blast_ResFreqStdComp(const BlastScoreBlk* sbp, Blast_ResFreq* rfp);

/** Fills in residue frequences for a given sequence.
 * @param sbp needed for alphabet information [in]
 * @param rfp object to be populated [in|out]
 * @param string sequence for calculation [in]
 * @param length length of above sequence [in]
 */
Int2
Blast_ResFreqString(const BlastScoreBlk* sbp, Blast_ResFreq* rfp, char* string, Int4 length);

Blast_ScoreFreq*
Blast_ScoreFreqNew(Int4 score_min, Int4 score_max);

Blast_ScoreFreq*
Blast_ScoreFreqFree(Blast_ScoreFreq* sfp);

/** Calculates the score frequencies.
 *
 * @param sbp object with scoring information [in]
 * @param sfp object to hold frequency information [in|out]
 * @param rfp1 letter frequencies for first sequence (query) [in]
 * @param rfp2 letter frequencies for second sequence (database) [in]
 * @return zero on success
 */
Int2
BlastScoreFreqCalc(const BlastScoreBlk* sbp, Blast_ScoreFreq* sfp,
		Blast_ResFreq* rfp1, Blast_ResFreq* rfp2);

/** Fill in the matrix for blastn using the penaly and rewards
 * The query sequence alphabet is blastna, the subject sequence
 * is ncbi2na.  The alphabet blastna is defined in blast_stat.h
 * and the first four elements of blastna are identical to ncbi2na.
 * if sbp->matrix==NULL, it is allocated.
 * @param sbp the BlastScoreBlk on which reward, penalty, and matrix will be set
 [in|out]
 * @return zero on success.
*/
Int2 BlastScoreBlkNuclMatrixCreate(BlastScoreBlk* sbp);

Int2
BlastScoreBlkMaxScoreSet(BlastScoreBlk* sbp);


#endif	/* BLAST_STAT_FUNCTIONS_H */

