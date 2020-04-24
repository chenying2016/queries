#include "search_setup.h"

#include "../../ncbi_blast/setup/blast_encoding.h"

/* See description in blast_setup.h */
Int2
Blast_ScoreBlkKbpGappedCalc(BlastScoreBlk * sbp,
                            const BlastScoringOptions * scoring_options,
                            EBlastProgramType program, 
                            const BlastQueryInfo * query_info,
                            Blast_Message** error_return)
{
    Int4 index = 0;
    Int2 retval = 0;

    if (sbp == NULL || scoring_options == NULL) {
        Blast_PerrorWithLocation(error_return, BLASTERR_INVALIDPARAM, -1);
        return 1;
    }

    /* Fill values for gumbel parameters*/
    if (program == eBlastTypeBlastn) {
        /* TODO gumbel parameters are not supported for nucleotide search yet  
        retval = 
                Blast_KarlinBlkNuclGappedCalc(sbp->kbp_gap_std[index],
                    scoring_options->gap_open, scoring_options->gap_extend, 
                    scoring_options->reward, scoring_options->penalty, 
                    sbp->kbp_std[index], &(sbp->round_down), error_return); */
    } else if (sbp->gbp) {
        retval = Blast_GumbelBlkCalc(sbp->gbp,
                    scoring_options->gap_open, scoring_options->gap_extend,
                    sbp->name, error_return);
    }
    if (retval)  return retval;

    /* Allocate and fill values for a gapped Karlin block, given the scoring
       options, then copy it for all the query contexts, as long as they're
       contexts that will be searched (i.e.: valid) */
    for (index = query_info->first_context;
         index <= query_info->last_context; index++) {

        if ( !query_info->contexts[index].is_valid ) {
            continue;
        }

        sbp->kbp_gap_std[index] = Blast_KarlinBlkNew();

        /* At this stage query sequences are nucleotide only for blastn */
        if (program == eBlastTypeBlastn) {
          /* If reward/penalty are both zero the calling program is
           * indicating that a matrix must be used to score both the
           * ungapped and gapped alignments.  If this is the case
           * set reward/penalty to allowed values so that extraneous
           * KA stats can be performed without error. -RMH-
           */
            if ( scoring_options->reward == 0 &&  scoring_options->penalty == 0 )
            {
              retval =
                Blast_KarlinBlkNuclGappedCalc(sbp->kbp_gap_std[index],
                    scoring_options->gap_open, scoring_options->gap_extend,
                    BLAST_REWARD, BLAST_PENALTY,
                    sbp->kbp_std[index], &(sbp->round_down), error_return);
            }else {
              retval =
                Blast_KarlinBlkNuclGappedCalc(sbp->kbp_gap_std[index],
                    scoring_options->gap_open, scoring_options->gap_extend,
                    scoring_options->reward, scoring_options->penalty,
                    sbp->kbp_std[index], &(sbp->round_down), error_return);
            }
        } else {
            retval = 
                Blast_KarlinBlkGappedCalc(sbp->kbp_gap_std[index],
                    scoring_options->gap_open, scoring_options->gap_extend,
                    sbp->name, error_return);
        }
        if (retval) {
            return retval;
        }

        /* For right now, copy the contents from kbp_gap_std to 
         * kbp_gap_psi (as in old code - BLASTSetUpSearchInternalByLoc) */
        if (program != eBlastTypeBlastn && program != eBlastTypeMapping) {
            sbp->kbp_gap_psi[index] = Blast_KarlinBlkNew();
            Blast_KarlinBlkCopy(sbp->kbp_gap_psi[index], 
                                sbp->kbp_gap_std[index]);
        }
    }

    /* Set gapped Blast_KarlinBlk* alias */
    sbp->kbp_gap = Blast_QueryIsPssm(program) ? 
        sbp->kbp_gap_psi : sbp->kbp_gap_std;

    return 0;
}

/** Fills a scoring block structure for a PHI BLAST search. 
 * @param sbp Scoring block structure [in] [out]
 * @param options Scoring options structure [in]
 * @param blast_message Structure for reporting errors [out]
 * @param get_path callback function for matrix path [in]
 */
static Int2
s_PHIScoreBlkFill(BlastScoreBlk* sbp, const BlastScoringOptions* options,
   Blast_Message** blast_message, GET_MATRIX_PATH get_path)
{
   Blast_KarlinBlk* kbp;
   char buffer[1024];
   Int2 status = 0;
   int index;

   sbp->read_in_matrix = TRUE;
   if ((status = Blast_ScoreBlkMatrixFill(sbp, get_path)) != 0)
      return status;
   kbp = sbp->kbp_gap_std[0] = Blast_KarlinBlkNew();
   /* Point both non-allocated Karlin block arrays to kbp_gap_std. */
   sbp->kbp_gap = sbp->kbp_gap_std;

   /* For PHI BLAST, the H value is not used, but it is not allowed to be 0, 
      so set it to 1. */
   kbp->H = 1.0;

   /* This is populated so that the checks for valid contexts don't fail,
    * note that this field is not used at all during a PHI-BLAST search */
   sbp->sfp[0] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);

   /* Ideal Karlin block is filled unconditionally. */
   status = Blast_ScoreBlkKbpIdealCalc(sbp);
   if (status)
      return status;

   if (0 == strcmp("BLOSUM62", options->matrix)) {
      kbp->paramC = 0.50;
      if ((11 == options->gap_open) && (1 == options->gap_extend)) {
         kbp->Lambda = 0.270;
         kbp->K = 0.047;
      } else if ((9 == options->gap_open) && (2 == options->gap_extend)) {
         kbp->Lambda = 0.285;
         kbp->K = 0.075;
      } else if ((8 == options->gap_open) && (2 == options->gap_extend)) {
         kbp->Lambda = 0.265;
         kbp->K = 0.046;
      } else if ((7 == options->gap_open) && (2 == options->gap_extend)) {
         kbp->Lambda = 0.243;
         kbp->K = 0.032;
      } else if ((12 == options->gap_open) && (1 == options->gap_extend)) {
         kbp->Lambda = 0.281;
         kbp->K = 0.057;
      } else if ((10 == options->gap_open) && (1 == options->gap_extend)) {
         kbp->Lambda = 0.250;
         kbp->K = 0.033;
      } else {
          status = -1;
      }
   } else if (0 == strcmp("PAM30", options->matrix)) { 
       kbp->paramC = 0.30;
       if ((9 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.295;
           kbp->K = 0.13;
       } else if ((7 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.306;
           kbp->K = 0.15;
       } else if ((6 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.292;
           kbp->K = 0.13;
       } else if ((5 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.263;
           kbp->K = 0.077;
       } else if ((10 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.309;
           kbp->K = 0.15;
       } else if ((8 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.270;
           kbp->K = 0.070;
       } else {
           status = -1;
       }
   } else if (0 == strcmp("PAM70", options->matrix)) { 
       kbp->paramC = 0.35;
       if ((10 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.291;
           kbp->K = 0.089;
       } else if ((8 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.303;
           kbp->K = 0.13;
       } else if ((7 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.287;
           kbp->K = 0.095;
       } else if ((6 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.269;
           kbp->K = 0.079;
       } else if ((11 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.307;
           kbp->K = 0.13;
       } else if ((9 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.269;
           kbp->K = 0.058;
       } else {
           status = -1;
       }
   } else if (0 == strcmp("BLOSUM80", options->matrix)) { 
       kbp->paramC = 0.40;
       if ((10 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.300;
           kbp->K = 0.072;
       } else if ((8 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.308;
           kbp->K = 0.089;
       } else if ((7 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.295;
           kbp->K = 0.077;
       } else if ((6 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.271;
           kbp->K = 0.051;
       } else if ((11 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.314;
           kbp->K = 0.096;
       } else if ((9 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.277;
           kbp->K = 0.046;
       } else {
           status = -1;
       }
   } else if (0 == strcmp("BLOSUM45", options->matrix)) { 
       kbp->paramC = 0.60;
       if ((14 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.199;
           kbp->K = 0.040;
       } else if ((13 == options->gap_open) && (3 == options->gap_extend)) {
           kbp->Lambda = 0.209;
           kbp->K = 0.057;
       } else if ((12 == options->gap_open) && (3 == options->gap_extend)) {
           kbp->Lambda = 0.203;
            kbp->K = 0.049;
       } else if ((11 == options->gap_open) && (3 == options->gap_extend)) {
           kbp->Lambda = 0.193;
           kbp->K = 0.037;
       } else if ((10 == options->gap_open) && (3 == options->gap_extend)) {
           kbp->Lambda = 0.182;
           kbp->K = 0.029;
       } else if ((15 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.206;
           kbp->K = 0.049;
       } else if ((13 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.190;
           kbp->K = 0.032;
       } else if ((12 == options->gap_open) && (2 == options->gap_extend)) {
           kbp->Lambda = 0.177;
           kbp->K = 0.023;
       } else if ((19 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.209;
           kbp->K = 0.049;
       } else if ((18 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.202;
           kbp->K = 0.041;
       } else if ((17 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.195;
           kbp->K = 0.034;
       } else if ((16 == options->gap_open) && (1 == options->gap_extend)) {
           kbp->Lambda = 0.183;
           kbp->K = 0.024;
       } else {
           status = -1;
       }
   } else {
       status = -2;
   }

   if (status == -1) {
       sprintf(buffer, "The combination %d for gap opening cost and %d for "
               "gap extension is not supported in PHI-BLAST with matrix %s\n",
               options->gap_open, options->gap_extend, options->matrix);
   } else if (status == -2) {
       sprintf(buffer, "Matrix %s not allowed in PHI-BLAST\n", options->matrix);
   }
   if (status) 
       Blast_MessageWrite(blast_message, eBlastSevWarning, kBlastMessageNoContext, buffer);
   else {

       /* fill in the rest of kbp_gap_std */
       for(index=1;index<sbp->number_of_contexts;index++)
       sbp->kbp_gap_std[index] = (Blast_KarlinBlk*)
           BlastMemDup(sbp->kbp_gap_std[0], sizeof(Blast_KarlinBlk));

       /* copy kbp_gap_std to kbp_std */
       for(index=0;index<sbp->number_of_contexts;index++)
       sbp->kbp_std[index] = (Blast_KarlinBlk*)
           BlastMemDup(sbp->kbp_gap_std[0], sizeof(Blast_KarlinBlk));

       sbp->kbp = sbp->kbp_std;
   }

   return status;
}

Int2
Blast_ScoreBlkMatrixInit(EBlastProgramType program_number, 
                  const BlastScoringOptions* scoring_options,
                  BlastScoreBlk* sbp,
                  GET_MATRIX_PATH get_path)
{
    Int2 status = 0;

    if ( !sbp || !scoring_options ) {
        return 1;
    }

    /* Matrix only scoring is used to disable the greedy extension 
       optimisations which avoid use of a full-matrix.  This is 
       currently only turned on in RMBlastN -RMH-  */
    sbp->matrix_only_scoring = FALSE;

    if (program_number == eBlastTypeBlastn ||
        program_number == eBlastTypeMapping) {

        BLAST_ScoreSetAmbigRes(sbp, 'N');
        BLAST_ScoreSetAmbigRes(sbp, '-');

        /* If reward/penalty are both zero the calling program is
         * indicating that a matrix must be used to score both the
         * ungapped and gapped alignments.  Set the new 
         * matrix_only_scoring.  For now reset reward/penalty to 
         * allowed blastn values so that extraneous KA stats can be 
         * performed without error. -RMH-
         */
        if ( scoring_options->penalty == 0 && scoring_options->reward == 0 )
        {
           sbp->matrix_only_scoring = TRUE;
           sbp->penalty = BLAST_PENALTY;
           sbp->reward = BLAST_REWARD;
        }else {
           sbp->penalty = scoring_options->penalty;
           sbp->reward = scoring_options->reward;
        }

        if (scoring_options->matrix && *scoring_options->matrix != NULLB) {
 
            sbp->read_in_matrix = TRUE;
            sbp->name = strdup(scoring_options->matrix);
 
        } else {
            char buffer[50];
            sbp->read_in_matrix = FALSE;
            sprintf(buffer, "blastn matrix:%ld %ld",
                    (long) sbp->reward, (long) sbp->penalty);
            sbp->name = strdup(buffer);
        }
 
     } else {
        sbp->read_in_matrix = TRUE;
        BLAST_ScoreSetAmbigRes(sbp, 'X');
        sbp->name = BLAST_StrToUpper(scoring_options->matrix);
    }
    status = Blast_ScoreBlkMatrixFill(sbp, get_path);
    if (status) {
        return status;
    }

    return status;
}

static Int2
s_JumperScoreBlkFill(BlastScoreBlk* sbp, const BlastQueryInfo* query_info,
                     Blast_Message** error_return)
{
    Int4 context;
    Blast_KarlinBlk* kbp;
    Int2 status;

    /* Create ungapped block */
    status = Blast_ScoreBlkKbpIdealCalc(sbp);
    if (status) {
        return status;
    }

    for (context = query_info->first_context;
         context <= query_info->last_context; ++context) {

        if (!query_info->contexts[context].is_valid) {
            continue;
        }

        sbp->sfp[context] = NULL;
        sbp->kbp_std[context] = Blast_KarlinBlkNew();
        Blast_KarlinBlkCopy(sbp->kbp_std[context], sbp->kbp_ideal);
    }
    sbp->kbp = sbp->kbp_std;

    /* Create gapped block */
    context = query_info->first_context;
    while (!query_info->contexts[context].is_valid) {
        context++;
    }

    sbp->kbp_gap_std[context] = kbp = Blast_KarlinBlkNew();
    status = Blast_KarlinBlkNuclGappedCalc(kbp, BLAST_GAP_OPEN_MEGABLAST,
                                           BLAST_GAP_EXTN_MEGABLAST,
                                           BLAST_REWARD,
                                           BLAST_PENALTY,
                                           sbp->kbp_std[context],
                                           &(sbp->round_down),
                                           error_return);
    if (status){
        return status;
    }

    for (++context;context <= query_info->last_context; ++context) {

        if (!query_info->contexts[context].is_valid) {
            continue;
        }

        sbp->kbp_gap_std[context] = Blast_KarlinBlkNew();
        Blast_KarlinBlkCopy(sbp->kbp_gap_std[context], kbp);
    }
    sbp->kbp_gap = sbp->kbp_gap_std;

    return status;
}

Int2 
BlastSetup_ScoreBlkInit(BLAST_SequenceBlk* query_blk, 
                        const BlastQueryInfo* query_info, 
                        const BlastScoringOptions* scoring_options, 
                        EBlastProgramType program_number, 
                        BlastScoreBlk* *sbpp, 
                        double scale_factor, 
                        Blast_Message* *blast_message,
                        GET_MATRIX_PATH get_path)
{
    BlastScoreBlk* sbp;
    Int2 status=0;      /* return value. */
    ASSERT(blast_message);

    if (sbpp == NULL)
       return 1;

    if (program_number == eBlastTypeBlastn ||
        program_number == eBlastTypeMapping) {

       sbp = BlastScoreBlkNew(BLASTNA_SEQ_CODE, query_info->last_context + 1);
       /* disable new FSC rules for nucleotide case for now */
       if (sbp && sbp->gbp) {
           sfree(sbp->gbp);
           sbp->gbp = NULL;
       }
    } else {
       sbp = BlastScoreBlkNew(BLASTAA_SEQ_CODE, query_info->last_context + 1);
    }

    if (!sbp) {
       Blast_PerrorWithLocation(blast_message, BLASTERR_MEMORY, -1);
       return 1;
    }

    *sbpp = sbp;
    sbp->scale_factor = scale_factor;

    /* Flag to indicate if we are using cross_match-like complexity
       adjustments on the raw scores.  RMBlastN is the currently 
       the only program using this flag. -RMH- */
    sbp->complexity_adjusted_scoring = scoring_options->complexity_adjusted_scoring;

    status = Blast_ScoreBlkMatrixInit(program_number, scoring_options, sbp, get_path);
    if (status) {
        Blast_PerrorWithLocation(blast_message, status, kBlastMessageNoContext);
        return status;
    }

    /* Fills in block for gapped blast. */
    if (Blast_ProgramIsPhiBlast(program_number)) {
       status = s_PHIScoreBlkFill(sbp, scoring_options, blast_message, get_path);
    } else if (Blast_ProgramIsMapping(program_number)) {
        /* Create fake score blocks for each query without computing base
           frequencies. We do not compute e-values for mapping, so the
           KA statistics are not needed, but the data structures are checked
           and used in BLAST engine. */
        status = s_JumperScoreBlkFill(sbp, query_info, blast_message);
    } else {
       status = Blast_ScoreBlkKbpUngappedCalc(program_number, sbp, query_blk->sequence, 
               query_info, blast_message);

       if (scoring_options->gapped_calculation) {
          status = 
              Blast_ScoreBlkKbpGappedCalc(sbp, scoring_options, program_number,
                                          query_info, blast_message);
       } else {
          ASSERT(sbp->kbp_gap == NULL);
          /* for ungapped cases we do not have gbp filled */
          if (sbp->gbp) {
              sfree(sbp->gbp);
              sbp->gbp=NULL;
          }
       }
    }

    return status;
}

Int2
BlastSetup_Validate(const BlastQueryInfo* query_info, 
                    const BlastScoreBlk* score_blk) 
{
    int index;
    Boolean valid_context_found = FALSE;
    ASSERT(query_info);

    for (index = query_info->first_context;
         index <= query_info->last_context;
         index++) {
        if (query_info->contexts[index].is_valid) {
            valid_context_found = TRUE;
        } else if (score_blk) {
            ASSERT(score_blk->kbp[index] == NULL);
            ASSERT(score_blk->sfp[index] == NULL);
            if (score_blk->kbp_gap) {
                ASSERT(score_blk->kbp_gap[index] == NULL);
            }
        }
    }

    if (valid_context_found) {
        return 0;
    } else {
        return 1;
    }
}

Int2 BLAST_MainSetUp(EBlastProgramType program_number,
    const QuerySetUpOptions *qsup_options,
    const BlastScoringOptions *scoring_options,
    BLAST_SequenceBlk *query_blk,
    const BlastQueryInfo *query_info,
    double scale_factor,
    BlastSeqLoc **lookup_segments, 
    BlastMaskLoc **mask,
    BlastScoreBlk **sbpp, 
    Blast_Message **blast_message,
    GET_MATRIX_PATH get_path)
{
    Int2 status = 0;            /* return value */
    status = BlastSetup_ScoreBlkInit(query_blk, query_info, scoring_options, 
                                     program_number, sbpp, scale_factor, 
                                     blast_message, get_path);
    hbn_assert(status == 0);
    if (status) {
        return status;
    }

    if ( (status = BlastSetup_Validate(query_info, *sbpp) != 0)) {
        if (*blast_message == NULL) {
            Blast_PerrorWithLocation(blast_message, status, kBlastMessageNoContext);
        }
        return 1;
    }

    return status;
}

BlastScoreBlk*
CSetupFactory__CreateScoreBlock(const HbnOptionsHandle* opts_memento,
                                BLAST_SequenceBlk* queries,
                                BlastQueryInfo* query_info)
{
    hbn_assert(opts_memento);

    double rps_scale_factor = 1.0;
    Blast_Message* blast_msg = NULL;

    BlastScoreBlk* retval = NULL;
    Int2 status = BLAST_MainSetUp(eBlastTypeBlastn,
                                  opts_memento->m_QueryOpts,
                                  opts_memento->m_ScoringOpts,
                                  queries,
                                  query_info,
                                  rps_scale_factor,
                                  NULL,
                                  NULL,
                                  &retval,
                                  &blast_msg,
                                  NULL);
    hbn_assert(status == 0);

    return retval;
}

/** Return the search space appropriate for a given context from
 * a list of tabulated search spaces
 * @param eff_len_options Container for search spaces [in]
 * @param context_index Identifier for the search space to return [in]
 * @param blast_message List of messages, to receive possible warnings [in][out]
 * @return The selected search space
 */
static Int8 s_GetEffectiveSearchSpaceForContext(
                        const BlastEffectiveLengthsOptions* eff_len_options,
                        int context_index, Blast_Message **blast_message)
{
    Int8 retval = 0;

    if (eff_len_options->num_searchspaces == 0) {
        retval = 0;
    } else if (eff_len_options->num_searchspaces == 1) {
        if (context_index != 0) {
            Blast_MessageWrite(blast_message, eBlastSevWarning, context_index, 
                    "One search space is being used for multiple sequences");
        }
        retval = eff_len_options->searchsp_eff[0];
    } else if (eff_len_options->num_searchspaces > 1) {
        ASSERT(context_index < eff_len_options->num_searchspaces);
        retval = eff_len_options->searchsp_eff[context_index];
    } else {
        abort();    /* should never happen */
    }
    return retval;
}

Int2 BLAST_CalcEffLengths (EBlastProgramType program_number, 
   const BlastScoringOptions* scoring_options,
   const BlastEffectiveLengthsParameters* eff_len_params, 
   const BlastScoreBlk* sbp, BlastQueryInfo* query_info,
   Blast_Message * *blast_message)
{
   double alpha=0, beta=0; /*alpha and beta for new scoring system */
   Int4 index;		/* loop index. */
   Int4	db_num_seqs;	/* number of sequences in database. */
   Int8	db_length;	/* total length of database. */
   Blast_KarlinBlk* *kbp_ptr; /* Array of Karlin block pointers */
   const BlastEffectiveLengthsOptions* eff_len_options = eff_len_params->options;

   if (!query_info || !sbp)
      return -1;


   /* use overriding value from effective lengths options or the real value
      from effective lengths parameters. */
   if (eff_len_options->db_length > 0)
      db_length = eff_len_options->db_length;
   else
      db_length = eff_len_params->real_db_length;

   /* If database (subject) length is not available at this stage, and
    * overriding value of effective search space is not provided by user,
    * do nothing.
    * This situation can occur in the initial set up for a non-database search,
    * where each subject is treated as an individual database. 
    */
   if (db_length == 0 &&
       !BlastEffectiveLengthsOptions_IsSearchSpaceSet(eff_len_options)) {
      return 0;
   }

   if (Blast_SubjectIsTranslated(program_number))
      db_length = db_length/3;  

   if (eff_len_options->dbseq_num > 0)
      db_num_seqs = eff_len_options->dbseq_num;
   else
      db_num_seqs = eff_len_params->real_num_seqs;
   
   /* Mapping does not need length correction */
   if (Blast_ProgramIsMapping(program_number)) {
        for (index = query_info->first_context;
           index <= query_info->last_context;
           index++) {
           query_info->contexts[index].eff_searchsp = db_length;
        }

        return 0;
   }

   /* PHI BLAST search space calculation is different. */
   if (Blast_ProgramIsPhiBlast(program_number))
   {
        for (index = query_info->first_context;
           index <= query_info->last_context;
           index++) {
           Int8 effective_search_space = db_length - (db_num_seqs*(query_info->contexts[index].length_adjustment));
           query_info->contexts[index].eff_searchsp = effective_search_space;
        }

        return 0;
   }

   /* N.B.: the old code used kbp_gap_std instead of the kbp_gap alias (which
    * could be kbp_gap_psi), hence we duplicate that behavior here */
   kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp);
   
   for (index = query_info->first_context;
        index <= query_info->last_context;
        index++) {
      Blast_KarlinBlk *kbp; /* statistical parameters for the current context */
      Int4 length_adjustment = 0; /* length adjustment for current iteration. */
      Int4 query_length;   /* length of an individual query sequence */
      
      /* Effective search space for a given sequence/strand/frame */
      Int8 effective_search_space =
          s_GetEffectiveSearchSpaceForContext(eff_len_options, index,
                                              blast_message);

      kbp = kbp_ptr[index];
      
      if (query_info->contexts[index].is_valid &&
          ((query_length = query_info->contexts[index].query_length) > 0) ) {

         /* Use the correct Karlin block. For blastn, two identical Karlin
          * blocks are allocated for each sequence (one per strand), but we
          * only need one of them.
          */
         if (program_number == eBlastTypeBlastn) {
             /* Setting reward and penalty to zero is being used to indicate
              * that matrix scoring should be used for ungapped and gapped
              * alignment.  For now reward/penalty are being reset to the
              * default blastn values to not disturb the KA calcs  -RMH- */
             if ( scoring_options->reward == 0 && scoring_options->penalty == 0 )
             {
                 Blast_GetNuclAlphaBeta(BLAST_REWARD,
                                    BLAST_PENALTY,
                                    scoring_options->gap_open,
                                    scoring_options->gap_extend,
                                    sbp->kbp_std[index],
                                    scoring_options->gapped_calculation,
                                    &alpha, &beta);
             }else {
                 Blast_GetNuclAlphaBeta(scoring_options->reward,
                                    scoring_options->penalty,
                                    scoring_options->gap_open,
                                    scoring_options->gap_extend,
                                    sbp->kbp_std[index],
                                    scoring_options->gapped_calculation,
                                    &alpha, &beta);
             }
         } else {
             BLAST_GetAlphaBeta(sbp->name, &alpha, &beta,
                                scoring_options->gapped_calculation, 
                                scoring_options->gap_open, 
                                scoring_options->gap_extend, 
                                sbp->kbp_std[index]);
         }
         BLAST_ComputeLengthAdjustment(kbp->K, kbp->logK,
                                       alpha/kbp->Lambda, beta,
                                       query_length, db_length,
                                       db_num_seqs, &length_adjustment);

         if (effective_search_space == 0) {

             /* if the database length was specified, do not
                adjust it when calculating the search space;
                it's counter-intuitive to specify a value and
                not have that value be used */
        	 /* Changing this rule for now sicne cutoff score depends
        	  * on the effective seach space length. SB-902
        	  */

        	 Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);

        	 // Just in case effective_db_length < 0
        	 if (effective_db_length <= 0)
        		 effective_db_length = 1;

             effective_search_space = effective_db_length *
                             (query_length - length_adjustment);
         }
      }
      query_info->contexts[index].eff_searchsp = effective_search_space;
      query_info->contexts[index].length_adjustment = length_adjustment;
   }
   return 0;
}

Int2 
BLAST_GapAlignSetUp(EBlastProgramType program_number,
    const text_t* db,
    const BlastScoringOptions* scoring_options,
    const BlastEffectiveLengthsOptions* eff_len_options,
    const BlastExtensionOptions* ext_options,
    const BlastHitSavingOptions* hit_options,
    const BlastInitialWordOptions* word_options,
    BlastQueryInfo* query_info, 
    BlastScoreBlk* sbp, 
    BlastScoringParameters** score_params,
    BlastExtensionParameters** ext_params,
    BlastHitSavingParameters** hit_params,
    BlastEffectiveLengthsParameters** eff_len_params,
    BlastInitialWordParameters** word_params)
{
   Int2 status = 0;
    Int8 total_length = seqdb_size(db);
    Int4 num_seqs = seqdb_num_seqs(db);
    Int8 avg_subject_length = total_length / num_seqs;

   /* Initialize the effective length parameters with real values of
      database length and number of sequences */
   BlastEffectiveLengthsParametersNew(eff_len_options, total_length, num_seqs, 
                                      eff_len_params);
   /* Effective lengths are calculated for all programs except PHI BLAST. */
   if ((status = BLAST_CalcEffLengths(program_number, scoring_options, 
                     *eff_len_params, sbp, query_info, NULL)) != 0)
   {
      *eff_len_params = BlastEffectiveLengthsParametersFree(*eff_len_params);
      return status;
   }

   if((status=BlastScoringParametersNew(scoring_options, sbp, score_params)) != 0)
   {
      *eff_len_params = BlastEffectiveLengthsParametersFree(*eff_len_params);
      *score_params = BlastScoringParametersFree(*score_params); 
      return status;
   }

   if((status=BlastExtensionParametersNew(program_number, ext_options, sbp, 
                               query_info, ext_params)) != 0)
   {
      *eff_len_params = BlastEffectiveLengthsParametersFree(*eff_len_params);
      *score_params = BlastScoringParametersFree(*score_params); 
      *ext_params = BlastExtensionParametersFree(*ext_params); 
      return status;
   }

   BlastHitSavingParametersNew(program_number, hit_options, sbp, query_info, 
                               avg_subject_length, (*ext_params)->options->compositionBasedStats, 
				hit_params);

    BlastInitialWordParametersNew(program_number, word_options,
        *hit_params, NULL, sbp, query_info,
        avg_subject_length, word_params);
   return status;
}