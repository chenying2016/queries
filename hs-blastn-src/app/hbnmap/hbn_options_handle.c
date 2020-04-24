#include "hbn_options_handle.h"

static void
HbnOptionsHandle_SetDustFiltering(HbnOptionsHandle* opts, Boolean val)
{
    if (opts->m_QueryOpts->filtering_options->dustOptions) {
        opts->m_QueryOpts->filtering_options->dustOptions = 
            SDustOptionsFree(opts->m_QueryOpts->filtering_options->dustOptions);
    }
    if (val == FALSE) return;
    SDustOptionsNew(&opts->m_QueryOpts->filtering_options->dustOptions);
}

static void
HbnOptionsHandle_SetQueryOptionDefaults(HbnOptionsHandle* opts)
{
    HbnOptionsHandle_SetDustFiltering(opts, TRUE);
    opts->m_QueryOpts->filtering_options->mask_at_hash = TRUE;
    opts->m_QueryOpts->strand_option = F_R;
}

static void
HbnOptionsHandle_SetLookupTableDefaults(HbnOptionsHandle* opts)
{
    opts->m_LutOpts->lut_type = eNaLookupTable;
    opts->m_LutOpts->word_size = BLAST_WORDSIZE_NUCL;
    opts->m_LutOpts->threshold = BLAST_WORD_THRESHOLD_BLASTN;
    opts->m_LutOpts->stride = 0;
}

static void
HbnOptionsHandle_SetMBLookupTableDefaults(HbnOptionsHandle* opts)
{
    opts->m_LutOpts->lut_type = eMBLookupTable;
    opts->m_LutOpts->word_size = BLAST_WORDSIZE_MEGABLAST;
    opts->m_LutOpts->threshold = BLAST_WORD_THRESHOLD_MEGABLAST;
    opts->m_LutOpts->stride = 0;
}

static void
HbnOptionsHandle_SetInitialWordOptionsDefaults(HbnOptionsHandle* opts)
{
    opts->m_InitWordOpts->x_dropoff = BLAST_UNGAPPED_X_DROPOFF_NUCL;
    opts->m_InitWordOpts->window_size = BLAST_WINDOW_SIZE_NUCL;
    opts->m_InitWordOpts->scan_range = BLAST_SCAN_RANGE_NUCL;
}

static void
HbnOptionsHandle_SetMBInitialWordOptionsDefaults(HbnOptionsHandle* opts)
{
    opts->m_InitWordOpts->window_size = BLAST_WINDOW_SIZE_NUCL;
}

static void
HbnOptionsHandle_SetGappedExtensionDefaults(HbnOptionsHandle* opts)
{
    opts->m_ExtnOpts->gap_x_dropoff = BLAST_GAP_X_DROPOFF_NUCL;
    opts->m_ExtnOpts->gap_x_dropoff_final = BLAST_GAP_X_DROPOFF_FINAL_NUCL;
    opts->m_InitWordOpts->gap_trigger = BLAST_GAP_TRIGGER_NUCL;
    opts->m_ExtnOpts->ePrelimGapExt = eDynProgScoreOnly;
    opts->m_ExtnOpts->eTbackExt = eDynProgTbck;
}

static void
HbnOptionsHandle_SetMBGappedExtensionDefaults(HbnOptionsHandle* opts)
{
    opts->m_ExtnOpts->gap_x_dropoff = BLAST_GAP_X_DROPOFF_GREEDY;
    opts->m_ExtnOpts->gap_x_dropoff_final = BLAST_GAP_X_DROPOFF_FINAL_NUCL;
    opts->m_InitWordOpts->gap_trigger = BLAST_GAP_TRIGGER_NUCL;
    opts->m_ExtnOpts->ePrelimGapExt = eGreedyScoreOnly;
    opts->m_ExtnOpts->eTbackExt = eGreedyTbck;
}

static void
HbnOptionsHandle_SetScoringOptionsDefaults(HbnOptionsHandle* opts)
{
    opts->m_ScoringOpts->matrix = NULL;
    opts->m_ScoringOpts->gap_open = BLAST_GAP_OPEN_NUCL;
    opts->m_ScoringOpts->gap_extend = BLAST_GAP_EXTN_NUCL;
    opts->m_ScoringOpts->reward = 2;
    opts->m_ScoringOpts->penalty = -3;
    opts->m_ScoringOpts->gapped_calculation = 1;
    opts->m_ScoringOpts->complexity_adjusted_scoring = 0;
    opts->m_ScoringOpts->is_ooframe = 0;
    opts->m_ScoringOpts->shift_pen = INT2_MAX;
}

static void
HbnOptionsHandle_SetMBScoringOptionsDefaults(HbnOptionsHandle* opts)
{
    opts->m_ScoringOpts->matrix = NULL;
    opts->m_ScoringOpts->gap_open = BLAST_GAP_OPEN_MEGABLAST;
    opts->m_ScoringOpts->gap_extend = BLAST_GAP_EXTN_MEGABLAST;
    opts->m_ScoringOpts->reward = 1;
    opts->m_ScoringOpts->penalty = -2;
    opts->m_ScoringOpts->gapped_calculation = 1;
    opts->m_ScoringOpts->complexity_adjusted_scoring = 0;
    opts->m_ScoringOpts->is_ooframe = 0;
    opts->m_ScoringOpts->shift_pen = INT2_MAX;
}

static void
HbnOptionsHandle_SetHitSavingOptionsDefaults(HbnOptionsHandle* opts)
{
    opts->m_HitSaveOpts->hitlist_size = 500;
    opts->m_HitSaveOpts->expect_value = BLAST_EXPECT_VALUE;
    opts->m_HitSaveOpts->percent_identity = 0;
    opts->m_HitSaveOpts->hsp_num_max = 0;
    opts->m_HitSaveOpts->max_hsps_per_subject = 0;
    opts->m_HitSaveOpts->min_diag_separation = 50;
    opts->m_HitSaveOpts->mask_level = 101;
    opts->m_HitSaveOpts->cutoff_score = 0;
    opts->m_HitSaveOpts->low_score_perc = 0;
    opts->m_HitSaveOpts->query_cov_hsp_perc = 0;
}

static void
HbnOptionsHandle_SetMBHitSavingOptionsDefaults(HbnOptionsHandle* opts)
{
    opts->m_HitSaveOpts->hitlist_size = 500;
    opts->m_HitSaveOpts->expect_value = BLAST_EXPECT_VALUE;
    opts->m_HitSaveOpts->percent_identity = 0;
    opts->m_HitSaveOpts->hsp_num_max = 0;
    opts->m_HitSaveOpts->max_hsps_per_subject = 0;
    opts->m_HitSaveOpts->min_diag_separation = 6;
    opts->m_HitSaveOpts->mask_level = 101;
    opts->m_HitSaveOpts->cutoff_score = 0;
    opts->m_HitSaveOpts->low_score_perc = 0;
    opts->m_HitSaveOpts->query_cov_hsp_perc = 0;
}

void
HbnOptionsHandle_SetEffectiveSearchSpace(HbnOptionsHandle* opts, Int8 eff)
{
    if (opts->m_EffLenOpts->num_searchspaces < 1) {
        opts->m_EffLenOpts->num_searchspaces = 1;
        if (opts->m_EffLenOpts->searchsp_eff) sfree(opts->m_EffLenOpts->searchsp_eff);
        opts->m_EffLenOpts->searchsp_eff = (Int8*)malloc(sizeof(Int8));
    }
    for (int i = 0; i < opts->m_EffLenOpts->num_searchspaces; ++i) {
        opts->m_EffLenOpts->searchsp_eff[i] = eff;
    }
}

static void
HbnOptionsHandle_SetEffectiveLengthsOptionsDefaults(HbnOptionsHandle* opts)
{
    opts->m_EffLenOpts->db_length = 0;
    opts->m_EffLenOpts->dbseq_num = 0;
    HbnOptionsHandle_SetEffectiveSearchSpace(opts, 0);
}

static void
HbnOptionsHandle_SetTraditionalBlastnDefaults(HbnOptionsHandle* opts)
{
    opts->m_Program = eBlastn;
    HbnOptionsHandle_SetQueryOptionDefaults(opts);
    HbnOptionsHandle_SetLookupTableDefaults(opts);
    HbnOptionsHandle_SetInitialWordOptionsDefaults(opts);
    HbnOptionsHandle_SetGappedExtensionDefaults(opts);
    HbnOptionsHandle_SetScoringOptionsDefaults(opts);
    HbnOptionsHandle_SetHitSavingOptionsDefaults(opts);
    HbnOptionsHandle_SetEffectiveLengthsOptionsDefaults(opts);
}

static void
HbnOptionsHandle_SetTraditionalMegablastDefaults(HbnOptionsHandle* opts)
{
    opts->m_Program = eMegablast;
    HbnOptionsHandle_SetQueryOptionDefaults(opts);
    HbnOptionsHandle_SetMBLookupTableDefaults(opts);
    HbnOptionsHandle_SetMBInitialWordOptionsDefaults(opts);
    HbnOptionsHandle_SetMBGappedExtensionDefaults(opts);
    HbnOptionsHandle_SetMBScoringOptionsDefaults(opts);
    HbnOptionsHandle_SetMBHitSavingOptionsDefaults(opts);
    HbnOptionsHandle_SetEffectiveLengthsOptionsDefaults(opts);
}

HbnOptionsHandle*
HbnOptionsHandleNew(EProgram program)
{
    HbnOptionsHandle* opts = (HbnOptionsHandle*)calloc(1, sizeof(HbnOptionsHandle));
    EBlastProgramType blast_program = eBlastTypeBlastn;
    BLAST_InitDefaultOptions(blast_program,
        &opts->m_LutOpts,
        &opts->m_QueryOpts,
        &opts->m_InitWordOpts,
        &opts->m_ExtnOpts,
        &opts->m_HitSaveOpts,
        &opts->m_ScoringOpts,
        &opts->m_EffLenOpts,
        &opts->m_PSIBlastOpts,
        &opts->m_DbOpts);
    if (program == eBlastn) HbnOptionsHandle_SetTraditionalBlastnDefaults(opts);
    else HbnOptionsHandle_SetTraditionalMegablastDefaults(opts);
    return opts;
}

HbnOptionsHandle*
HbnOptionsHandleFree(HbnOptionsHandle* opts)
{
    LookupTableOptionsFree(opts->m_LutOpts);
    BlastQuerySetUpOptionsFree(opts->m_QueryOpts);
    BlastInitialWordOptionsFree(opts->m_InitWordOpts);
    BlastExtensionOptionsFree(opts->m_ExtnOpts);
    BlastHitSavingOptionsFree(opts->m_HitSaveOpts);
    BlastScoringOptionsFree(opts->m_ScoringOpts);
    BlastEffectiveLengthsOptionsFree(opts->m_EffLenOpts);
    PSIBlastOptionsFree(opts->m_PSIBlastOpts);
    BlastDatabaseOptionsFree(opts->m_DbOpts);
    free(opts);
    return NULL;
}

EBlastProgramType
EProgramToEBlastProgramType(EProgram p)
{
    switch (p) {
    case eBlastn:
    case eMegablast:
    case eDiscMegablast:
    case eVecScreen:
        return eBlastTypeBlastn;

    case eMapper:
        return eBlastTypeMapping;
        
    case eBlastp:
        return eBlastTypeBlastp;
        
    case eBlastx:
        return eBlastTypeBlastx;
        
    case eTblastn:
        return eBlastTypeTblastn;
        
    case eTblastx:
        return eBlastTypeTblastx;
        
    case eRPSBlast:
        return eBlastTypeRpsBlast;
        
    case eRPSTblastn:
        return eBlastTypeRpsTblastn;
        
    case ePSIBlast:
    case eDeltaBlast:
        return eBlastTypePsiBlast;
        
    case ePSITblastn:
        return eBlastTypePsiTblastn;

    case ePHIBlastp:
        return eBlastTypePhiBlastp;
        
    case ePHIBlastn:
        return eBlastTypePhiBlastn;
        
    default:
        return eBlastTypeUndefined;
    }
}

void
HbnOptionsHandle_Update(const HbnProgramOptions* opts, HbnOptionsHandle* handle)
{
    /// general search options
    handle->m_Program = opts->task;
    handle->m_HitSaveOpts->expect_value = opts->expect_value;
    handle->m_ScoringOpts->reward = opts->reward;
    handle->m_ScoringOpts->penalty = opts->penalty;
    handle->m_ScoringOpts->gap_open = opts->gap_open;
    handle->m_ScoringOpts->gap_extend = opts->gap_extend;

    /// input query options
    if (opts->strand == FWD) {
        handle->m_QueryOpts->strand_option = 1;
    } else if (opts->strand == REV) {
        handle->m_QueryOpts->strand_option = 2;
    } else {
        hbn_assert(opts->strand == F_R);
        handle->m_QueryOpts->strand_option = 3;
    }

    /// query filtering options
    if (!opts->use_dust_masker) {
        if (handle->m_QueryOpts->filtering_options->dustOptions) {
            SDustOptionsFree(handle->m_QueryOpts->filtering_options->dustOptions);
            handle->m_QueryOpts->filtering_options->dustOptions = NULL;
        }
    } else {
        if (handle->m_QueryOpts->filtering_options->dustOptions) {
            SDustOptionsFree(handle->m_QueryOpts->filtering_options->dustOptions);
            handle->m_QueryOpts->filtering_options->dustOptions = NULL;
        }
        SDustOptionsNew(&handle->m_QueryOpts->filtering_options->dustOptions);
        handle->m_QueryOpts->filtering_options->dustOptions->level = opts->dust_level;
        handle->m_QueryOpts->filtering_options->dustOptions->window = opts->dust_window;
        handle->m_QueryOpts->filtering_options->dustOptions->linker = opts->dust_linker;        
    }

    /// restrict search or results
    handle->m_HitSaveOpts->percent_identity = opts->perc_identity;
    handle->m_HitSaveOpts->query_cov_hsp_perc = opts->query_cov_hsp_perc;
    handle->m_HitSaveOpts->max_hsps_per_subject = opts->max_hsps_per_subject;
    handle->m_HitSaveOpts->hitlist_size = opts->hitlist_size;

    if (opts->keep_best_hsp_per_subject) {
        if (handle->m_HitSaveOpts->hsp_filt_opt == NULL) {
            handle->m_HitSaveOpts->hsp_filt_opt = BlastHSPFilteringOptionsNew();
        }
        if (handle->m_HitSaveOpts->hsp_filt_opt->subject_besthit_opts == NULL) {
            BOOL isProtein = Blast_ProgramIsNucleotide(EProgramToEBlastProgramType(handle->m_Program));
            BlastHSPSubjectBestHitOptions* besthit = BlastHSPSubjectBestHitOptionsNew(isProtein);
            BlastHSPFilteringOptions_AddSubjectBestHit(handle->m_HitSaveOpts->hsp_filt_opt, &besthit);
            hbn_assert(besthit == NULL);
        }
    } 

    /// statistical options
    if (opts->searchsp_eff) HbnOptionsHandle_SetEffectiveSearchSpace(handle, opts->searchsp_eff);
    if (opts->db_length) handle->m_EffLenOpts->db_length = opts->db_length;
    if (opts->dbseq_num) handle->m_EffLenOpts->dbseq_num = opts->dbseq_num;

    Blast_Message* blast_msg = NULL;
    if (BLAST_ValidateOptions(eBlastTypeBlastn,
        handle->m_ExtnOpts,
        handle->m_ScoringOpts,
        handle->m_LutOpts,
        handle->m_InitWordOpts,
        handle->m_HitSaveOpts,
        &blast_msg)) {
        HBN_ERR("Options validate fail");
    }
}

HbnProgramOptions*
HbnProgramOptionsFree(HbnProgramOptions* opts)
{
    if (opts->db_dir) free((void*)opts->db_dir);
    opts->db_dir = NULL;

    if (opts->output) free((void*)opts->output);
    opts->output = NULL;

    free(opts);

    return NULL;
}