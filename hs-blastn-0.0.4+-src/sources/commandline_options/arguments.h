#ifndef BLAST_ARGS_H
#define	BLAST_ARGS_H

#include "options.h"

#include <vector>

void PrintHelpSimple(bool note_line);
void PrintHelpFull();

struct SingleCmdLineArg
{
    int (*ArgDealFunctionPtr)(void* dst, int& i, int argc, const char** argv);
    void* arg_ptr;
    const char* arg_name;
};

struct CmdLineArgs
{
    PrintVersionOption* version_options;
    HelpOptions* help_options;
    SBlastFilterOptions* filtering_options;
    BlastInitialWordOptions* word_options;
    BlastExtensionOptions* ext_options;
    BlastHitSavingOptions* hit_options;
    BlastEffectiveLengthsOptions* eff_options;
    BlastScoringOptions* scoring_options;
    SeedingOptions* seed_options;
    InputOptions* input_options;
    OutputOptions* output_options;
    RunningOptions* running_options;

    std::vector<SingleCmdLineArg> cmd_args;

    CmdLineArgs(PrintVersionOption* po,
                HelpOptions* ho, 
                SBlastFilterOptions* fo,
                BlastInitialWordOptions* iwo,
                BlastExtensionOptions* eo,
                BlastHitSavingOptions* hso,
                BlastEffectiveLengthsOptions* elo,
                BlastScoringOptions* so,
                SeedingOptions* seedo,
                InputOptions* io,
                OutputOptions* oo,
                RunningOptions* ro)
            :   version_options(po), 
                help_options(ho),
                filtering_options(fo),
                word_options(iwo),
                ext_options(eo),
                hit_options(hso),
                eff_options(elo),
                scoring_options(so),
                seed_options(seedo),
                input_options(io),
                output_options(oo),
                running_options(ro)
    {}

    ~CmdLineArgs()
    {}

    void InitialArgList();
    void ParseCmdLineArgs(int argc, const char** argv);
    void PrintArgsNames();

    void PrintSimpleUsage();
    void PrintFullUsage();
};

struct Options
{
    SBlastFilterOptions* filtering_options;
    BlastInitialWordOptions* word_options;
    BlastExtensionOptions* ext_options;
    BlastHitSavingOptions* hit_options;
    //BlastEffectiveLengthsOptions* eff_options;
    BlastScoringOptions* scoring_options;
    SeedingOptions* seed_options;
    InputOptions* input_options;
    OutputOptions* output_options;
    RunningOptions* running_options;
    PrintVersionOption* version_options;
    HelpOptions* help_options;
    
    CmdLineArgs* cmd_args;

    Options();
    ~Options();
    void ParseCmdLineArgs(int argc, const char** argv);
};

#endif	/* BLAST_ARGS_H */

