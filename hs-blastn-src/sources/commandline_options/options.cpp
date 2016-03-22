#include "options.h"
#include "stat_functions.h"
#include "index.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
using std::clog;

PrintVersionOption*
PrintVersionOptionNew()
{
    PrintVersionOption* ret = 
        (PrintVersionOption*)calloc(sizeof( PrintVersionOption), 1);
    ret->print_version = FALSE;

    return ret;
}

HelpOptions*
HelpOptionsNew()
{
    HelpOptions* ret = 
        (HelpOptions*)calloc(sizeof(HelpOptions), 1);
    ret->simple_help = FALSE;
    ret->full_help = FALSE;

    return ret;
}

SeedingOptions*
SeedingOptionsNew()
{
    SeedingOptions* seed_options = 
            (SeedingOptions*)calloc(sizeof(SeedingOptions), 1);
    seed_options->seed_size = 28;
    
    return seed_options;
}

InputOptions*
InputOptionsNew()
{
    InputOptions* options = 
            (InputOptions*)calloc(sizeof(InputOptions), 1);
    options->db = NULL;
    options->query = NULL;
	options->query_list = NULL;
    
    return options;
}

OutputOptions*
OutputOptionsNew()
{
    OutputOptions* options = 
            (OutputOptions*)calloc(sizeof(OutputOptions), 1);
    options->num_alignments = 250;
    options->num_descriptions = 500;
    options->out = stdout;
    options->output_file_name = NULL;
    options->outfmt = ePairwise;
    return options;
}

RunningOptions*
RunningOptionsNew()
{
    RunningOptions* options = 
            (RunningOptions*)calloc(sizeof(RunningOptions), 1);
    options->num_threads = 1;
    options->query_strand = eStrandBoth;
    return options;
}

BlastInitialWordOptions*
BlastInitialWordOptionsNew()
{
	BlastInitialWordOptions* options =
			(BlastInitialWordOptions*)calloc(sizeof(BlastInitialWordOptions), 1);
	options->gap_trigger = 27.0;
	options->window_size = 0;
	options->scan_range = 0;
	options->x_dropoff = 0.0;

	return options;
}

BlastExtensionOptions*
BlastExtensionOptionsNew()
{
	BlastExtensionOptions* options =
			(BlastExtensionOptions*)calloc(sizeof(BlastExtensionOptions), 1);
	options->gap_x_dropoff = 25.0;
	options->gap_x_dropoff_final = 100.0;
	options->ePrelimGapExt = eGreedyScoreOnly;
	options->compositionBasedStats = FALSE;

	return options;
}

BlastHitSavingOptions*
BlastHitSavingOptionsNew()
{
	BlastHitSavingOptions* options =
			(BlastHitSavingOptions*)calloc(sizeof(BlastHitSavingOptions), 1);
	options->expect_value = 10.0;
	options->cutoff_score = 0;
	options->percent_identity = 0.0;
	options->hitlist_size = 500;
	options->hsp_num_max = 0;
	options->total_hsp_limit = 0;
	options->culling_limit = 0;
	options->mask_level = 101;
	options->min_hit_length = 0;
	options->min_diag_separation = 6;

	return options;
}

BlastEffectiveLengthsOptions*
BlastEffectiveLengthsOptionsNew()
{
	BlastEffectiveLengthsOptions* options =
			(BlastEffectiveLengthsOptions*)calloc(sizeof(BlastEffectiveLengthsOptions), 1);
	options->db_length = 0;
	options->dbseq_num = 0;
	options->num_searchspaces = 1;
	//options->searchsp_eff = (Int8*)calloc(sizeof(Int8), 1);
        options->searchsp_eff = NULL;

	return options;
}

BlastScoringOptions*
BlastScoringOptionsNew()
{
	BlastScoringOptions* options =
			(BlastScoringOptions*)calloc(sizeof(BlastScoringOptions), 1);

	options->gap_open = 0;
	options->gap_extend = 0;
	options->reward = 1;
	options->penalty = -2;
	options->gapped_calculation =  TRUE;

	return options;
}

void PrintScoringOptions(BlastScoringOptions* options)
{
	cerr << "[" << __func__ << "] reward = " << options->reward << endl;
	cerr << "[" << __func__ << "] penalty = " << options->penalty << endl;
	cerr << "[" << __func__ << "] gap open = " << options->gap_open << endl;
	cerr << "[" << __func__ << "] gap extend = " << options->gap_extend << endl;
}
void
PrintInitialWordOptions(BlastInitialWordOptions* options)
{
	cerr << "[" << __func__ << "] gap_trigger = " << options->gap_trigger << endl;
	cerr << "[" << __func__ << "] window_size = " << options->window_size << endl;
	cerr << "[" << __func__ << "] scan_range = " << options->scan_range << endl;
	cerr << "[" << __func__ << "] x_dropoff = " << options->x_dropoff << endl;
	cerr << endl << endl;
}

void
PrintExtensionOptions(BlastExtensionOptions* options)
{
	cerr << "[" << __func__ << "] gap_x_dropoff = " << options->gap_x_dropoff << endl;
	cerr << "[" << __func__ << "] gap_x_dropoff_final = " << options->gap_x_dropoff_final << endl;
	cerr << "[" << __func__ << "] ePrelimGapExt = " << options->ePrelimGapExt << endl;
	cerr << "[" << __func__ << "] compositionBasedStats = " << options->compositionBasedStats << endl;
	cerr << endl << endl;
}

void
PrintHitSavingOptions(BlastHitSavingOptions* options)
{
    fprintf(stderr, "[%s] expect_value = %f\n", __func__, options->expect_value);
    fprintf(stderr, "[%s] cutoff_score = %d\n", __func__, options->cutoff_score);
    fprintf(stderr, "[%s] percent_identity = %f\n", __func__, options->percent_identity);
    fprintf(stderr, "[%s] hitlist_size = %d\n", __func__, options->hitlist_size);
    fprintf(stderr, "[%s] hsp_num_max = %d\n", __func__, options->hsp_num_max);
    fprintf(stderr, "[%s] total_hsp_limit = %d\n", __func__, options->total_hsp_limit);
    fprintf(stderr, "[%s] culling_limit = %d\n", __func__, options->culling_limit);
    fprintf(stderr, "[%s] mask_level = %d\n", __func__, options->mask_level);
    fprintf(stderr, "[%s] min_hit_length = %d\n", __func__, options->min_hit_length);
    fprintf(stderr, "[%s] min_diag_separation = %d\n", __func__, options->min_diag_separation);
    fprintf(stderr, "\n\n");
}

void PrintEffectiveLengthOptions(BlastEffectiveLengthsOptions* options)
{
	fprintf(stderr, "[%s] db_length = %ld\n", __func__, (long)options->db_length);
	fprintf(stderr, "[%s] dbseq_num = %d\n", __func__, options->dbseq_num);
	fprintf(stderr, "[%s] num_searchspaces = %ld\n", __func__, (long)options->num_searchspaces);
	if (options->num_searchspaces > 0 && options->searchsp_eff != NULL)
	{
		int i;
		for (i = 0; i < options->num_searchspaces; ++i)
			fprintf(stderr, "[%s] searchsp_eff[%d] = %ld\n", __func__, i, (long)options->searchsp_eff[i]);
	}
	fprintf(stderr, "\n\n");
}

Int2 BlastExtensionOptionsValidate(const BlastExtensionOptions* options)
{
    if (options == NULL)
    {
        fprintf(stderr, "[%s] Error: The BlastExtensionOptions pointer cannot be Null.\n",
                __func__);
        exit(1);
    }
    return 1;
}

Int2 BlastScoringOptionsValidate(const BlastScoringOptions* options)
{
    if (options == NULL)
    {
        fprintf(stderr, "[%s] Error: The BlastScoringOptions pointer cannot be Null.\n",
                __func__);
        exit(1);
    }   
    if (!(options->penalty == 0 && options->reward == 0))
    {
        if (options->penalty >= 0)
        {
                fprintf(stderr, "[%s] Error: BLASTN penalty must be negative.\n",
                        __func__);
                exit(1);            
        }
        if (options->gapped_calculation && 
            !BLAST_CheckRewardPenaltyScores(options->reward, options->penalty))
        {
                fprintf(stderr, "[%s] Error: BLASTN reward/penalty combination not supported for gapped search.\n",
                        __func__);
                exit(1);              
        }
        if (options->gapped_calculation && options->gap_open > 0 && options->gap_extend == 0)
        {
                fprintf(stderr, "[%s] Error: BLASTN gap extension penalty cannot be 0.\n",
                        __func__);
                exit(1);            
        }
    }
    
    return 1;
}

Int2 BlastInitialWordOptionsValidate(const BlastInitialWordOptions* options)
{
    if (options == NULL)
    {
        fprintf(stderr, "[%s] Error: The BlastInitialWordOptions pointer cannot be Null.\n",
                __func__);
        exit(1);
    }      
    
    return 1;
}

Int2 BlastHitSavingOptionsValidate(const BlastHitSavingOptions* options)
{
    if (options == NULL)
    {
        fprintf(stderr, "[%s] Error: The BlastHitSavingOptions pointer cannot be Null.\n",
                __func__);
        exit(1);
    } 
    if (options->hitlist_size < 1)
    {
        fprintf(stderr, "[%s] Error: No hits are being saved.\n",
                __func__);
        exit(1);        
    }
    if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
    {
        fprintf(stderr, "[%s] Error: expect value or cutoff score must be greater than zero.\n",
                __func__);
        exit(1);        
    }
	if (options->percent_identity < 0.0 || options->percent_identity > 100.0)
	{
		fprintf(stderr, "[%s] Error: Argument \"perc_identity\". Illegal value, expected [0,100]: \'%g\'\n", __func__, options->percent_identity);
		exit(1);
	}
    
    return 1;
}

Int2 BlastExtensionScoringOptionsValidate(const BlastExtensionOptions* ext_options,
                                          const BlastScoringOptions* score_options)
{
    if (ext_options == NULL)
    {
        fprintf(stderr, "[%s] Error: The BlastExtensionOptions pointer cannot be Null.\n",
                __func__);
        exit(1);
    } 
    if (score_options == NULL)
    {
        fprintf(stderr, "[%s] Error: The BlastScoringOptions pointer cannot be Null.\n",
                __func__);
        exit(1);
    }   
    
    return 1;
}

Int2 SeedingOptionsValidate(const SeedingOptions* options)
{
    if (options == NULL)
    {
        fprintf(stderr, "[%s] Error: The SeedingOptions pointer cannot be Null.\n",
                __func__);
        exit(1);
    }     
    if (options->seed_size < FMIndex::kLutSize)
    {
        fprintf(stderr, "[%s] Error: Word-size must be %d or greater.\n",
                __func__, FMIndex::kLutSize);
        exit(1);
    }    
    
    return 1;
}

Int2 InputOutputOptionsValidate(const InputOptions* ioptions,
                                const OutputOptions* ooptions)
{
    if (ioptions->db == NULL)
    {
        fprintf(stderr, "[%s] Error: No Database is specified.\n",
                __func__);
        exit(1);        
    }
    if (ioptions->query == NULL && ioptions->query_list == NULL)
    {
        fprintf(stderr, "[%s] Error: No Query is specified.\n",
                __func__);
        exit(1);          
    }
    if ((ooptions->outfmt != 0) &&
       (ooptions->outfmt != 6) &&
       (ooptions->outfmt != 7))
    {
        fprintf(stderr, "[%s] Error: Currently only formats 0, 6, 7 are supported.\n",
                __func__);
        exit(1);          
    }        
    
    return 1;
}

Int2 RunningOptionsValidate(const RunningOptions* options)
{
    if (options == NULL)
    {
        fprintf(stderr, "[%s] Error: The RunningOptions pointer cannot be Null.\n",
                __func__);
        exit(1);
    }  
    if (options->num_threads < 1)
    {
        fprintf(stderr, "[%s] Error: Number of threads must be positive.\n",
                __func__);
        exit(1);        
    }
	
	int physical_cpus = sysconf(_SC_NPROCESSORS_CONF);
	if (options->num_threads > physical_cpus)
	{
		fprintf(stderr, "\n\n[%s] Warning: You specify %d searching threads, but only %d CPU cores are dectected.\n\n",
				__func__, options->num_threads, physical_cpus);
	}
	
    return 1;
}

Int2 VersionOptionsValidate(const PrintVersionOption* version_options)
{
    if (version_options->print_version)
    {
        fprintf(stderr, "hs-blastn: version 0.0.5+\n");
        return 1;
    }

    return 0;
}

Int2 HelpOptionsValidate(const HelpOptions* help_options)
{
    extern void PrintHelpSimple(bool);
    extern void PrintHelpFull();

    if (help_options->simple_help)
    {
        PrintHelpSimple(true);
        return 1;
    }
    if (help_options->full_help)
    {
        PrintHelpSimple(false);
        PrintHelpFull();
        return 1;
    }

    return 0;
}

Int2 FilterOptionsValidate(const SBlastFilterOptions* filtering_options)
{
//    if (filtering_options->windowMaskerOptions->database == NULL)
//    {
//        fprintf(stderr, "[%s] Error: The -window_masker_db option must be specified.\n", __func__);
//        exit(1);
//   }

    return 1;
}

Int2 BLAST_ValidateOptions(const BlastExtensionOptions* ext_options,
                           const BlastScoringOptions* score_options,
                           const BlastInitialWordOptions* word_options,
                           const BlastHitSavingOptions* hit_options,
                           const InputOptions* ioptions,
                           const OutputOptions* ooptions,
                           const SeedingOptions* seed_options,
                           const RunningOptions* running_options,
                           const SBlastFilterOptions* filtering_options)
{
    Int2 status = 0;

    InputOutputOptionsValidate(ioptions, ooptions);
    
    BlastExtensionOptionsValidate(ext_options);
    
    BlastScoringOptionsValidate(score_options);
    
    BlastInitialWordOptionsValidate(word_options);
    
    BlastHitSavingOptionsValidate(hit_options);
    
    BlastExtensionScoringOptionsValidate(ext_options, score_options);
    
    SeedingOptionsValidate(seed_options);
    
    RunningOptionsValidate(running_options);

    FilterOptionsValidate(filtering_options);
    
    return status;
}
