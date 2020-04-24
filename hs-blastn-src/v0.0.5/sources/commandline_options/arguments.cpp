#include "arguments.h"
#include "utility.h"

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <iostream>
#include <sstream>

/// option names

/// input query options
static const char* kQuery               = "-query";
static const char* kStrand              = "-strand";
static const char* kQueryList           = "-query_list";

/// general search options
static const char* kDb                  = "-db";
static const char* kOut                 = "-out";
static const char* kEvalue              = "-evalue";
static const char* kWordSize            = "-word_size";
static const char* kGapOpen             = "-gapopen";
static const char* kGapExtend           = "-gapextend";
static const char* kPenalty				= "-penalty";
static const char* kReward				= "-reward";

/// blast-2-sequences options
static const char* kSubject				= "-subject";

/// formatting options
static const char* kOutFmt				= "-outfmt";
static const char* kNumDescriptions		= "-num_descriptions";
static const char* kNumAlignments		= "-num_alignments";
static const char* kMaxTargetSequences	= "-max_target_seqs";

/// query filtering options
static const char* kDust                = "-dust";
static const char* kWinowMaskerDb       = "-window_masker_db";

/// extension options
static const char* kXdropUngap          = "-xdrop_ungap";
static const char* kXdropGap            = "-xdrop_gap";
static const char* kXdropGapFinal       = "-xdrop_gap_final";
static const char* kMinRawGappedScore   = "-min_raw_gapped_score";
static const char* kPercentIdentity     = "-perc_identity";

/// miscellaneous options
static const char* kNumThreads          = "-num_threads";
static const char* kVersion             = "-version";
static const char* kSimpleHelp          = "-h";
static const char* kFullHelp            = "-help";

/// end options names

static bool
ValidateIntArguments(const char* src)
{
    int len = strlen(src);
    int i;

    if (len == 0) return false;

    if (src[0] == '-')
    {
        if (len == 1) return false;
        for (i = 1; i < len; ++i)
        if (!isdigit(src[i]))
            return false;
    }
    else
    {
        for (i = 0; i < len; ++i)
        if (!isdigit(src[i]))
            return false;
    }

    return true;
}

static int Int4ArgDealFunction(void* dst, int &i, int argc, const char** argv)
{
    if (i >= argc) 
    {
        fprintf(stderr, "Error: No value is specified for %s\n", argv[i-1]);
        exit(1);
    }
    if (!ValidateIntArguments(argv[i])) return -2;

    Int4* arg_ptr = (Int4*)dst;
    *arg_ptr = atoi(argv[i]);
    ++i;
    return 1;
}

static int Int2ArgDealFunction(void* dst, int &i, int argc, const char** argv)
{
    if (i >= argc) 
    {
        fprintf(stderr, "Error: No value is specified for %s\n", argv[i-1]);
        exit(1);
    }
    if (!ValidateIntArguments(argv[i])) return -2;

    Int2* arg_ptr = (Int2*)dst;
    *arg_ptr = atoi(argv[i]);
    ++i;
    return 1;   
}

static double StringToDouble(const char* str)
{
    if (str == NULL)
    {
		cy_utility::Log::ErrorAndExit(NULL, "Error: No double value is supplied.");
    }
    char c = str[0];
    if (!isdigit((unsigned int)c) && (c != '.') && (c != '-') && (c != '+'))
    {
		cy_utility::Log::ErrorAndExit(NULL, "Error: %s is not a Double number.", str);
    }

    // conversion
    char* endptr = 0;
    const char* begptr = str;
    double n;

    n = strtod(begptr, &endptr);

    return n;
}

static int FloatArgDealFunction(void* dst, int& i, int argc, const char** argv)
{
    double* ptr = (double*)dst;
    if (i >= argc)
    {
        fprintf(stderr, "Error, No value is specified for %s.\n", argv[i-1]);
        exit(1);
    }
    *ptr = StringToDouble(argv[i]);
    ++i;
    return 1;
}

static int BooleanArgDealFunction(void* dst, int& i, int argc, const char** argv)
{
    Boolean* ptr = (Boolean*)dst;
    *ptr = TRUE;
    return 1;
}

static int DustArgDealFunction(void* dst, int& i, int argc, const char** argv)
{
    SBlastFilterOptions* fo = (SBlastFilterOptions*)dst;

	if (i >= argc)
	{
		cy_utility::Log::ErrorAndExit(cy_utility::kAligner, "no argument(s) is supplied to option -dust");
	}
    if (strcmp(argv[i], "yes") == 0)
    {
        fo->mask_at_seeding = TRUE;
		++i;
        return 1;
    }
    else if (strcmp(argv[i], "no") == 0)
    {
        fo->mask_at_seeding = FALSE;
		++i;
        return 1;
    }
    else
    {
        if (argc - i < 3)
        {
            fprintf(stderr, "Error, Arguments to Dust options are not correct.\n");
            exit(1);
        }
        int n;
        if (!ValidateIntArguments(argv[i]))
        {
            fprintf(stderr, "Error, Invalid number of arguments to option");
            fprintf(stderr, " %s\n", kDust);
            exit(1);
        }
        n = atoi(argv[i]);
        fo->dustOptions->level = n;
        if (!ValidateIntArguments(argv[i + 1]))
        {
            fprintf(stderr, "Error, Invalid number of arguments to option");
            fprintf(stderr, " %s\n", kDust);
            exit(1);
        }
        n = atoi(argv[i + 1]);
        fo->dustOptions->linker = n;
        if (!ValidateIntArguments(argv[i + 2]))
        {
            fprintf(stderr, "Error, Invalid number of arguments to option");
            fprintf(stderr, " %s\n", kDust);
            exit(1);
        }
        n = atoi(argv[i + 2]);
        fo->dustOptions->window = n;
        i += 3;
        return 1;
    }
}

static int StringArgDealFunction(void* dst, int& i, int argc, const char** argv)
{
    if (i >= argc)
    {
        fprintf(stderr, "Error in parsing arguments: Please specified a value.\n");
        exit(1);
    }    
    
    const char** str = (const char**)dst;
    *str = argv[i];
    ++i;

    return 1;
}

static int StrandArgDealFunction(void* dst, int& i, int argc, const char** argv)
{
    ESequenceStrand* strand = (ESequenceStrand*)dst;

    if (i >= argc)
    {
        fprintf(stderr, "Error in parsing arguments: Please specified a strand.\n");
        exit(1);
    }

    const char* src_strand = argv[i];

    if (strcmp(src_strand, "plus") == 0)
        *strand = eStrandPlus;
    else if (strcmp(src_strand, "minus") == 0)
        *strand = eStrandMinus;
    else if (strcmp(src_strand, "both") == 0)
        *strand = eStrandBoth;
    else
    {
        fprintf(stderr, "Error in parsing arguments: %s is not valid for argument %s",
                kStrand, src_strand);
        exit(1);
    }
    ++i;

    return 1;
}

static void FillSingleCmeLineArg(SingleCmdLineArg& arg,
                          int (*funptr)(void* dst, int& i, int argc, const char** argv),
                          const char* arg_name,
                          void* arg_ptr)
{
    arg.ArgDealFunctionPtr = funptr;
    arg.arg_name = arg_name;
    arg.arg_ptr = arg_ptr;
}

std::vector<SingleCmdLineArg>::iterator
CmdLineArgsFindName(std::vector<SingleCmdLineArg>& cmd_args, const char* name)
{
    std::vector<SingleCmdLineArg>::iterator iter;
    for (iter = cmd_args.begin(); iter != cmd_args.end(); ++iter)
        if (strcmp(iter->arg_name, name) == 0)
            break;
    return iter;
}

void CmdLineArgs::InitialArgList()
{
    SingleCmdLineArg arg;
    cmd_args.clear();

    FillSingleCmeLineArg(arg, BooleanArgDealFunction, kVersion, &version_options->print_version);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, BooleanArgDealFunction, kSimpleHelp, &help_options->simple_help);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, BooleanArgDealFunction, kFullHelp, &help_options->full_help);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, StringArgDealFunction, kQuery, &input_options->query);
    cmd_args.push_back(arg);
	
	FillSingleCmeLineArg(arg, StringArgDealFunction, kQueryList, &input_options->query_list);
	cmd_args.push_back(arg);

    //FillSingleCmeLineArg(arg, StrandArgDealFunction, kStrand, &running_options->query_strand);
    //cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, StringArgDealFunction, kDb, &input_options->db);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, StringArgDealFunction, kOut, &output_options->output_file_name);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, FloatArgDealFunction, kEvalue, &hit_options->expect_value);///
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int4ArgDealFunction, kWordSize, &seed_options->seed_size);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int4ArgDealFunction, kGapOpen, &scoring_options->gap_open);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int4ArgDealFunction, kGapExtend, &scoring_options->gap_extend);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int2ArgDealFunction, kPenalty, &scoring_options->penalty);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int2ArgDealFunction, kReward, &scoring_options->reward);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int4ArgDealFunction, kOutFmt, &output_options->outfmt);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int4ArgDealFunction, kNumDescriptions, &output_options->num_descriptions);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int4ArgDealFunction, kNumAlignments, &output_options->num_alignments);
    cmd_args.push_back(arg);
	
	FillSingleCmeLineArg(arg, Int4ArgDealFunction, kMaxTargetSequences, &hit_options->hitlist_size);
	cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, DustArgDealFunction, kDust, filtering_options);
    cmd_args.push_back(arg);
    
    FillSingleCmeLineArg(arg, StringArgDealFunction, kWinowMaskerDb, 
                         &filtering_options->windowMaskerOptions->database);
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, FloatArgDealFunction, kXdropUngap, &word_options->x_dropoff); ///
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, FloatArgDealFunction, kXdropGap, &ext_options->gap_x_dropoff); ///
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, FloatArgDealFunction, kXdropGapFinal, &ext_options->gap_x_dropoff_final); ///
    cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int4ArgDealFunction, kMinRawGappedScore, &hit_options->cutoff_score);
    cmd_args.push_back(arg);
    
    FillSingleCmeLineArg(arg, FloatArgDealFunction, kPercentIdentity, &hit_options->percent_identity); ///
	cmd_args.push_back(arg);

    FillSingleCmeLineArg(arg, Int4ArgDealFunction, kNumThreads, &running_options->num_threads);
    cmd_args.push_back(arg);
}

void CmdLineArgs::PrintArgsNames()
{
    std::vector<SingleCmdLineArg>::iterator iter;
    for (iter = cmd_args.begin(); iter != cmd_args.end(); ++iter)
        fprintf(stderr, "%s\n", iter->arg_name);
}

void PrintHelpFull()
{
#define out std::cout
#define nline std::endl
static const char* kFourSpaceMargins = "    ";
static const char* kOneSPaceMargins = " ";
static const char* kSixSpaceMargins = "      ";

    out << "OPTIONAL ARGUMENTS" << nline;
    out << kOneSPaceMargins << "-h" << nline;
    out << kFourSpaceMargins << "Print USAGE and DESCRIPTION; ignore all other parameters" << nline;
    out << kOneSPaceMargins << "-help" << nline;
    out << kFourSpaceMargins << "Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters" << nline;
    out << kOneSPaceMargins << "-version" << nline;
    out << kFourSpaceMargins << "Print version number; ignore other arguments" << nline;
    out << nline;

    out << kOneSPaceMargins << "*** Input query options" << nline;
    out << kOneSPaceMargins;
    out << "-query <File_In>" << nline;
    out << kFourSpaceMargins;
    out << "Input file name (FASTA or FASTQ format)," << nline
		<< kFourSpaceMargins << "the base qualities in FASTQ formated files will be ignored." << nline;
	//out << "-query_list <File_In>" << nline;
	//out << kFourSpaceMargins;
	//out << "A file that contains the list of all the query files, one line for a file name." << nline;

    //out << kOneSPaceMargins;
    //out << "-strand <String, 'both', 'minus', 'plus'>" << nline;
    //out << kFourSpaceMargins;
    //out << "Query strand(s) to search against database" << nline;
    //out << kFourSpaceMargins;
    //out << "Default = 'both'" << nline;
    //out << nline;

    out << kOneSPaceMargins;
    out << "*** General search options" << nline;
    out << kOneSPaceMargins;
    out << "-db <File_In>" << nline;
    out << kFourSpaceMargins;
    out << "database name" << nline;
    out << kOneSPaceMargins;
    out << "-out <File_Out>" << nline;
    out << kFourSpaceMargins;
    out << "Output file name" << nline;
    out << kFourSpaceMargins;
    out << "Default = standard output" << nline;
    out << kOneSPaceMargins;
    out << "-evalue <Real>" << nline;
    out << kFourSpaceMargins;
    out << "Expectation value (E) threshold for saving hits" << nline;
    out << kFourSpaceMargins;
    out << "Default = '10'" << nline;
    out << kOneSPaceMargins;
    out << "-word_size <Integer, >= 12>" << nline;
    out << kFourSpaceMargins;
    out << "Word size for wordfiner algorithm (length of best perfect match)" << nline;
    out << kOneSPaceMargins;
    out << "-gapopen <Integer>" << nline;
    out << kFourSpaceMargins;
    out << "Cost to open a gap" << nline;
    out << kOneSPaceMargins;
    out << "-gapextend <Integer>" << nline;
    out << kFourSpaceMargins;
    out << "Cost to extend a gap" << nline;
    out << kOneSPaceMargins;
    out << "-penalty <Integer, <=0>" << nline;
    out << kFourSpaceMargins;
    out << "Penalty for a nucleotide mismatch" << nline;
    out << kOneSPaceMargins;
    out << "-reward <Interger, >=0>" << nline;
    out << kFourSpaceMargins;
    out << "Reward for a nucleotide match" << nline;
	out << "[Note] Not all the combinations [gapopen, gapextend, penalty, reward]" << nline 
		<< kFourSpaceMargins
		<< "   are supported by NCBI-BLASTN. And so is HS-BLASTN." << nline;
	out << nline;

    out << kOneSPaceMargins;
    out << "*** Formatting options" << nline;
    out << kOneSPaceMargins;
    out << "-outfmt <String>" << nline;
    out << kFourSpaceMargins;
    out << "alignment view options" << nline;
    out << kSixSpaceMargins;
    out << "0 = pairwise," << nline;
    out << kSixSpaceMargins;
    out << "6 = tabular," << nline;
    out << kSixSpaceMargins;
    out << "7 = tabular with comment lines" << nline;

    out << kOneSPaceMargins;
    out << "-num_descriptions <Integer, >=0>" << nline;
    out << kFourSpaceMargins;
    out << "Number of database sequences to show one-line descriptions for" << nline;
    out << kFourSpaceMargins;
    out << "Not applicable for outfmt >= 6" << nline;
    out << kFourSpaceMargins;
    out << "Default = '500'" << nline;
    out << kOneSPaceMargins;
    out << "-num_alignments <Integer, >=0>" << nline;
    out << kFourSpaceMargins;
    out << "Number of database sequences to show alignments for" << nline;
    out << kFourSpaceMargins;
    out << "Default = '250'" << nline;
    out << nline;

    out << kOneSPaceMargins;
    out << "-dust <String>" << nline;
    out << kFourSpaceMargins;
    out << "Filter query sequence with DUST (Format: 'yes', 'level window linker', or" << nline;
    out << kFourSpaceMargins;
    out << "'no' to disable)" << nline;
    out << kFourSpaceMargins;
    out << "Default = '20 64 1'" << nline;
    out << kOneSPaceMargins;
    out << "-perc_identity <Real, 0..100>" << nline;
    out << kFourSpaceMargins;
    out << "Percent identity" << nline;
    out << nline;
	out << "-max_target_seqs <Integer, >=1>" << nline;
	out << kFourSpaceMargins;
	out << "Maxinum number of aligned sequences to keep" << nline;
	out << kFourSpaceMargins;
	out << "Not applicable for outfmt <= 4" << nline;
	out << kFourSpaceMargins;
	out << "Default = '500'" << nline;
	out << kFourSpaceMargins;
	out << " * Incompatable with: num_descriptions, num_alignments" << nline;
    out << nline;

    out << kOneSPaceMargins;
    out << "*** Extension options" << nline;
    out << kOneSPaceMargins;
    out << "-xdrop_ungap <Real>" << nline;
    out << kFourSpaceMargins;
    out << "X-dropoff value (in bits) for ungapped extensions" << nline;
    out << kOneSPaceMargins;
    out << "-xdrop_gap <Real>" << nline;
    out << kFourSpaceMargins;
    out << "X-dropoff value (in bits) for preliminary gapped extensions" << nline;
    out << kOneSPaceMargins;
    out << "-xdrop_gap_final" << nline;
    out << kFourSpaceMargins;
    out << "X-dropoff value (in bits) for final gapped alignment" << nline;
    out << kOneSPaceMargins;
    out << "-min_raw_gapped_score <Integer>" << nline;
    out << kFourSpaceMargins;
    out << "Minimum raw gapped score to keep an alignment in the preliminary gapped and" << nline;
    out << kFourSpaceMargins;
    out << "traceback stages" << nline;
    out << nline;

    out << kOneSPaceMargins;
    out << "*** Miscellaneous options" << nline;
    out << kOneSPaceMargins;
    out << "-num_threads <Integer, >=1>" << nline;
    out << kFourSpaceMargins;
    out << "Number of threads (CPUs) to use in the search" << nline;
    out << kFourSpaceMargins;
    out << "Default = '1'" << nline;
}

void PrintHelpSimple(bool note_line)
{
static const char* kTwoSpaceMargins = "  ";
static const char* kFourSpaceMargins = "    ";
#define out std::cout
#define nline std::endl

   out << "USAGE" << nline;
   out << kTwoSpaceMargins;
   out << "hs-blastn align";
   out << " [-h] [-help] [-db database_name]" << nline;
   out << kFourSpaceMargins;
   out << "[-query input_file]" << nline;
   out << kFourSpaceMargins; 
   out << "[-evalue evalue] [-word_size int_value]" << nline;
   out << kFourSpaceMargins;
   out << "[-gapoptn open_penalty] [-gapextend extend_penalty]" << nline;
   out << kFourSpaceMargins;
   out << "[-perc_identity float_value] [-xdrop_ungap float_value]" << nline;
   out <<kFourSpaceMargins;
   out << "[-xdrop_gap float_value] [-xdrop_gap_final float_value]" << nline;
   out << kFourSpaceMargins;
   out << "[-penalty penalty] [-reward reward]" << nline;
   out << kFourSpaceMargins;
   out << "[-min_raw_gapped_score int_value] [-dust DUST_options]" << nline;
   //out << kFourSpaceMargins;
   //out << "[-strand strand] [-outfmt format]" << nline;
   out << kFourSpaceMargins;
   out << "[-num_descriptions int_value] [-num_alignments int_value]" << nline;
   out << kFourSpaceMargins;
   out << "[-num_threads int_value]" << nline;
   out << nline;
   out << kFourSpaceMargins;
   out << "[-max_target_seqs num_sequences]" << nline;
   out << nline;
   
    out << "DESCRIPTION" << nline;
    out << kFourSpaceMargins << "Nucleotide-Nucleotide Aligner" << nline;
    out << nline;
    
    if (note_line)
    {
        out << "Use '-help' to print detailed descriptions of command line arguments\n";
    }
}

void CmdLineArgs::ParseCmdLineArgs(int argc, const char** argv)
{
    static const char* simple_help = "-h";
    static const char* full_help = "-help";
    int i = 0;
    std::vector<SingleCmdLineArg>::iterator iter;
    while (i < argc)
    {
        //printf("%s %s\n", argv[i], argv[i+1]);
        iter = CmdLineArgsFindName(cmd_args, argv[i]);
        ++i;
        if (iter == cmd_args.end())
        {
            fprintf(stderr, "Fatal Error: Unrecognized arguments: %s\n",
                    argv[i]);
            exit(1);
        }
        int r = iter->ArgDealFunctionPtr(iter->arg_ptr, i, argc, argv);
        if (r != 1)
        {
            fprintf(stderr, "Fatal Error: Fail to parse %s\n",
                    argv[i - 1]);
            exit(1);
        }
    }
	
	bool num_desc_filled = false;
	bool num_align_filled = false;
	bool num_tseqs_filled = false;
	for (i = 0; i < argc; ++i)
	{
		if (strcmp(argv[i], kNumDescriptions) == 0) num_desc_filled = true;
		else if (strcmp(argv[i], kNumAlignments) == 0) num_align_filled = true;
		else if (strcmp(argv[i], kMaxTargetSequences) == 0) num_tseqs_filled = true;
	}
	
	if (output_options->outfmt <= eFlatQueryAnchoredIdentities)
	{
		if (num_tseqs_filled)
		{
			std::ostringstream err;
			err << "The parameter -max_target_seqs is ignored for "
				<< "output formats, 0,1,2,3. Use -num_descriptions "
				<< "and -num_alignments to control output";
			error_and_exit(err.str());
		}
	}
	else
	{
		if (num_desc_filled)
		{
			std::ostringstream err;
			err << "The parameter -num_descriptions is ignored for "
				<< "output formats > 4 . Use -max_target_seqs "
				<< "to control output";
			error_and_exit(err.str());
		}
		
		if (!num_tseqs_filled && num_align_filled)
			hit_options->hitlist_size = output_options->num_alignments;
	}
}

Options::Options()
{
    version_options = PrintVersionOptionNew();
    help_options = HelpOptionsNew();
    seed_options = SeedingOptionsNew();
    input_options = InputOptionsNew();
    output_options = OutputOptionsNew();
    running_options = RunningOptionsNew();
    word_options = BlastInitialWordOptionsNew();
    ext_options = BlastExtensionOptionsNew();
    hit_options = BlastHitSavingOptionsNew();
    //eff_options = BlastEffectiveLengthsOptionsNew();
    scoring_options = BlastScoringOptionsNew();
    filtering_options = new SBlastFilterOptions();
    
    BlastEffectiveLengthsOptions* eff_options = BlastEffectiveLengthsOptionsNew();
    
    cmd_args = new CmdLineArgs(version_options, help_options, 
                               filtering_options, word_options,
                               ext_options,
                               hit_options, eff_options, scoring_options,
                               seed_options, input_options, output_options,
                               running_options);
    free(eff_options);
}

Options::~Options()
{
    if (version_options) free(version_options);
    if (help_options) free(help_options);
    if (seed_options) free(seed_options);
    if (input_options) free(input_options);
    if (output_options) free(output_options);
    if (running_options) free(running_options);
    if (filtering_options) delete filtering_options;
    if (word_options) free(word_options);
    if (ext_options) free(ext_options);
    if (hit_options) free(hit_options);
    //if (eff_options) free(eff_options);
    if (scoring_options) free(scoring_options);
    
    if (cmd_args)
        delete cmd_args;
}

void Options::ParseCmdLineArgs(int argc, const char** argv)
{
    cmd_args->InitialArgList();
    cmd_args->ParseCmdLineArgs(argc, argv);

    Int2 status;
    status = VersionOptionsValidate(version_options);
    if (status) exit(0);

    status = HelpOptionsValidate(help_options);
    if (status) exit(0);
    
    BLAST_ValidateOptions(ext_options,
                          scoring_options,
                          word_options,
                          hit_options,
                          input_options,
                          output_options,
                          seed_options,
                          running_options,
                          filtering_options);
}




