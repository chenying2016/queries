#include "cmdline_args.h"

#include "../../corelib/hbn_aux.h"
#include "../../ncbi_blast/cmdline_args/blast_args.hpp"
#include "../../ncbi_blast/setup/blast_types.hpp"
#include "../../ncbi_blast/cmdline_args/cmdline_flags.hpp"

#include <limits>
#include <sstream>

BEGIN_NCBI_SCOPE
USING_SCOPE(blast);
USING_SCOPE(align_format);

using namespace std;

/// general search options
extern const string blast::kTask;
EProgram kDfltTask = eMegablast;
extern const string blast::kArgEvalue;
const double kDfltEvalue = BLAST_EXPECT_VALUE;
extern const string blast::kArgMismatch;
const int kDfltMismatch = MISMATCH_PENALTY;
extern const string blast::kArgMatch;
const int kDfltMatch = MATCH_REWARD;
extern const string blast::kArgGapOpen;
const int kDfltGapOpen = GAP_OPEN;
extern const string blast::kArgGapExtend;
const int kDfltGapExtend = GAP_EXTEND;
const string kArgSkipOverhang("skip_overhang");
const bool kDfltSkipOverhang = false;

/// input query options
extern const string blast::kArgStrand;
extern const string blast::kDfltArgStrand;

/// database options
const string kArgDbDir("db_dir");
const string kDfltDbDir("hbndb");
const string kArgKeepDb("keep_db");
const bool kDfltKeepDb = false;
const string kArgMinQuerySize("min_query_size");
const int kDfltMinQuerySize = 0;
const string kArgMaxQueryVolSeqs("max_query_seqs");
const int kDfltMaxQueryVolSeqs = 2000000;
const string kArgMaxQueryVolRes("max_query_res");
const size_t kDfltMaxQueryVolRes = static_cast<size_t>(4000000000);
const string kArgMinSubjectSize("min_subject_size");
const int kDfltMinSubjectSize = 0;
const string kArgMaxSubjectVolSeqs("max_subject_vol_seqs");
const int kDfltMaxSubjectVolSeqs = numeric_limits<int>::max();
const string kArgMaxSubjectVolRes("max_subject_vol_res");
const size_t kDfltMaxSubjectVolRes = static_cast<size_t>(4000000000);

/// ddf scoring options
const string kArgKmerSize("kmer_size");
const int kDfltKmerSize = 15;
const string kArgKmerWindow("kmer_window");
const int kDfltKmerWindow = 10;
const string kArgMaxKmerOcc("max_kmer_occ");
const int kDfltMaxKmerOcc = 1000;
const string kArgBlockSize("block_size");
const int kDfltBlockSize = 2000;
const string kArgDDFScore("ddf_score");
const int kDfltDDFScore = 2;

/// mem scoring options
const string kArgMemScKmerSize("chain_kmer_size");
const int kDfltMemScKmerSize = 13;
const int kMemScKmerSizeMin = 8;
const int kMemScKmerSizeMax = 13;
const string kArgMemScKmerWindow("chain_kmer_window");
const int kDfltMemScKmerWindow = 5;
const string kArgMemScMemScore("chain_score");
const int kDfltMemScMemScore = 13;

/// formatting options
extern const string align_format::kArgOutputFormat;
//extern const int align_format::kDfltArgOutputFormat;
EOutputFormat kDfltOutfmt = eSAM;
const string kArgSamRgInfo("sam_rg");
const char* kDfltSamRgInfo = NULL;
const string kArgSamRgSample("sam_sample");
const char* kDfltSamRgSample = NULL;
const string kArgSamDumpMd("sam_md");
const bool kDfltSamDumpMd = false;

/// query filtering options
extern const string blast::kArgDustFiltering;
extern const string blast::kDfltArgDustFiltering;
extern const string blast::kArgUseLCaseMasking;
extern const bool blast::kDfltArgUseLCaseMasking;

/// restrict search or results
extern const string blast::kArgPercentIdentity;
const double kDfltPercentIdentity = 0.0;
extern const string blast::kArgQueryCovHspPerc;
const double kDfltQueryCovHspPerc = 0.0;
const string kArgQueryCovHspRes("qcov_hsp_res");
const int kDfltQueryCovHspRes = 0;
extern const string blast::kArgMaxHSPsPerSubject;
const int kDfltMaxHSPsPerSubject = 10;
extern const string blast::kArgMaxTargetSequences;
const int kDfltMaxTargetSequences = BLAST_HITLIST_SIZE;
extern const string blast::kArgSubjectBestHit;
const bool kDfltSubjectBestHit = false;

/// statistical options
extern const string blast::kArgNumThreads;
extern const int blast::kDfltNumThreads;
extern const string blast::kArgGrid;
extern const int blast::kDfltNodeId;
extern const int blast::kDfltNumNodes;
extern const string blast::kArgOutput;
const string kDfltOutput("-");


/// groups
static const char* kGroupGeneralSearchOptions = "General search options";
static const char* kGroupInputQuery = "Input query options";
static const char* kGroupDbOptions = "Database options";
static const char* kGroupDDFSc = "DDF scoring options";
static const char* kGroupMemSc = "Chaining scoring options";
static const char* kGroupFormat = "Formatting options";
static const char* kGroupQueryFiltering = "Query Filtering Options";
static const char* kGroupRestrictSearch = "Restrict search or results";
static const char* kGroupStatistics = "Statistical options";
static const char* kGroupExtension = "Extension options";
static const char* kGroupMiscellaneous = "Miscellaneous options";

class CommandLineArguments : public IBlastCmdLineArgs
{
public:
    CommandLineArguments(HbnProgramOptions* options);
    ~CommandLineArguments();

    virtual void SetArgumentDescriptions(CArgDescriptions& arg_desc);

    virtual void ExtractAlgorithmOptions(const CArgs& cmd_line_args, CBlastOptions& options);

private:
    HbnProgramOptions*              m_Options;
};

CommandLineArguments::CommandLineArguments(HbnProgramOptions* options)
{
    m_Options = options;
}

CommandLineArguments::~CommandLineArguments()
{

}

void CommandLineArguments::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    /// create the groups so that the ordering is established
    arg_desc.SetCurrentGroup(kGroupInputQuery);
    arg_desc.SetCurrentGroup(kGroupGeneralSearchOptions);
    arg_desc.SetCurrentGroup(kGroupDbOptions);
    arg_desc.SetCurrentGroup(kGroupDDFSc);
    arg_desc.SetCurrentGroup(kGroupMemSc);
    arg_desc.SetCurrentGroup(kGroupFormat);
    arg_desc.SetCurrentGroup(kGroupQueryFiltering);
    arg_desc.SetCurrentGroup(kGroupRestrictSearch);
    arg_desc.SetCurrentGroup(kGroupStatistics);
    arg_desc.SetCurrentGroup(kGroupExtension);
    arg_desc.SetCurrentGroup(kGroupMiscellaneous);

    /// general search options
    arg_desc.SetCurrentGroup(kGroupGeneralSearchOptions);

    string kDefaultTask = EProgramToTaskName(kDfltTask);
    arg_desc.AddDefaultKey(kTask, "task_name",
                "Task to execute",
                CArgDescriptions::eString,
                kDefaultTask);
    set<string> supported_tasks = GetTasks(eNuclNucl);
    arg_desc.SetConstraint(kTask, new CArgAllowStringSet(supported_tasks));

    arg_desc.AddDefaultKey(kArgEvalue, "evalue",
                "Expectation value (E) threshold for saving hits",
                CArgDescriptions::eDouble,
                NStr::DoubleToString(kDfltEvalue));

    arg_desc.AddOptionalKey(kArgMismatch, "penalty",
                "Penalty for a nucleotide mismatch",
                CArgDescriptions::eInteger);

    arg_desc.AddOptionalKey(kArgMatch, "reward",
                "Reward for a nucleotide match",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMatch, new CArgAllowValuesGreaterThanOrEqual(0));
    
    arg_desc.AddOptionalKey(kArgGapOpen, "open_penalty",
                "Cost to open a gap",
                CArgDescriptions::eInteger);
    
    arg_desc.AddOptionalKey(kArgGapExtend, "extend_penalty",
                "Cost to extend a gap",
                CArgDescriptions::eInteger);

    /// input query options
    arg_desc.SetCurrentGroup(kGroupInputQuery);

    arg_desc.AddDefaultKey(kArgStrand, "strand",
                "Query strand(s) to search against database/subject",
                CArgDescriptions::eString,
                kDfltArgStrand);
    arg_desc.SetConstraint(kArgStrand, &(* new CArgAllow_Strings, kDfltArgStrand, "plus", "minus"));

    /// database options
    arg_desc.SetCurrentGroup(kGroupDbOptions);

    arg_desc.AddOptionalKey(kArgDbDir, "directory",
                "Directory to store the query and subject database",
                CArgDescriptions::eString);

    arg_desc.AddFlag(kArgKeepDb, "Do not delete the database after search?", true);

    arg_desc.AddOptionalKey(kArgMinQuerySize, "int_value",
                "Skip query sequences shorter than this value",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMinQuerySize, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgMaxQueryVolSeqs, "int_value",
                "Maximum number of sequences in one volume of query database",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMinQuerySize, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMaxQueryVolRes, "volume_size",
                "Maximum number of residues in one volume of query database",
                CArgDescriptions::eDataSize,
                NStr::UInt8ToString_DataSize(kDfltMaxQueryVolRes));
    arg_desc.SetConstraint(kArgMaxQueryVolRes, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddOptionalKey(kArgMinSubjectSize, "int_value",
                "Skip subject sequences shorter than this value",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMinSubjectSize, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgMaxSubjectVolSeqs, "int_value",
                "Maximum number of sequences in one volume of subject database",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMaxSubjectVolSeqs, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMaxSubjectVolRes, "volume_size",
                "Maximum number of residues in one volume of subject database",
                CArgDescriptions::eDataSize,
                NStr::UInt8ToString_DataSize(kDfltMaxSubjectVolRes));
    arg_desc.SetConstraint(kArgMaxSubjectVolRes, CArgAllowValuesGreaterThanOrEqual(1));

    /// DDF scoring options
    arg_desc.SetCurrentGroup(kGroupDDFSc);

    arg_desc.AddOptionalKey(kArgKmerSize, "int_value",
                "Kmer size for DDF scoring algorithm (length of best perfect match)",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgKmerSize, CArgAllowValuesBetween(10, 32));

    arg_desc.AddOptionalKey(kArgKmerWindow, "int_value",
                "Kmer sampling window size in subject sequences",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgKmerWindow, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgBlockSize, "int_value",
                "Split subject database into consecutive blocks, each having this number of residues",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgBlockSize, CArgAllowValuesGreaterThanOrEqual(500));

    arg_desc.AddOptionalKey(kArgDDFScore, "int_value",
                "Minimum DDF score of a candidate",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgDDFScore, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMaxKmerOcc, "int_value",
                "Filter out kmers occur larger than this value",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMaxKmerOcc));
    arg_desc.SetConstraint(kArgMaxKmerOcc, CArgAllowValuesGreaterThanOrEqual(1));

    //// mem scoring options
    arg_desc.SetCurrentGroup(kGroupMemSc);

    arg_desc.AddDefaultKey(kArgMemScKmerSize, "int_value",
                "Length of perfect matched kmers that are to be extended to MEMs",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMemScKmerSize));
    arg_desc.SetConstraint(kArgMemScKmerSize, CArgAllowValuesBetween(kMemScKmerSizeMin, kMemScKmerSizeMax));

    arg_desc.AddDefaultKey(kArgMemScKmerWindow, "int_value",
                "Kmer sampling window",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMemScKmerWindow));
    arg_desc.SetConstraint(kArgMemScKmerWindow, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMemScMemScore, "int_value",
                "Minimum chaining score of a candidate",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMemScMemScore));
    arg_desc.SetConstraint(kArgMemScMemScore, CArgAllowValuesGreaterThanOrEqual(0));

    /// output format
    arg_desc.SetCurrentGroup(kGroupFormat);

    string OutputFormatDescription = string(
        "alignment view options:\n"
        "  0 = Pairwise,\n"
        "  6 = Tabular,\n"
        "  7 = Tabular with comment lines,\n"
        " 17 = Sequence Alignment/Map (SAM)"
    );
    arg_desc.AddDefaultKey(kArgOutputFormat, "format",
                OutputFormatDescription,
                CArgDescriptions::eString,
                NStr::IntToString(kDfltArgOutputFormat));

    arg_desc.AddOptionalKey(kArgSamRgInfo, "rg_header",
                "Read group header line such as '@RG\tID:xyz\tSM:abc', only apply for SAM output format",
                CArgDescriptions::eString);

    arg_desc.AddOptionalKey(kArgSamRgSample, "sample_name",
                "Sample name for all read groups. Only apply for SAM output format",
                CArgDescriptions::eString);

    arg_desc.AddFlag(kArgSamDumpMd, "Add MD tag to the SAM output", true);

    /// query filtering options
    arg_desc.SetCurrentGroup(kGroupQueryFiltering);

    arg_desc.AddDefaultKey(kArgDustFiltering, "DUST_options",
                "Filter query sequence with DUST " 
                "(Format: '" + kDfltArgApplyFiltering + "', " +
                "'level window linker', or '" + kDfltArgNoFiltering +
                "' to disable)",
                CArgDescriptions::eString,
                kDfltArgDustFiltering);

    arg_desc.AddFlag(kArgUseLCaseMasking, 
                "Use lower case filtering in query and subject sequence(s)?", true);

    /// restrict search or results
    arg_desc.SetCurrentGroup(kGroupRestrictSearch);

    arg_desc.AddOptionalKey(kArgPercentIdentity, "float_value",
                "Percent identity",
                CArgDescriptions::eDouble);
    arg_desc.SetConstraint(kArgPercentIdentity, CArgAllow_Doubles(0.0, 100.0));

    arg_desc.AddOptionalKey(kArgQueryCovHspPerc, "float_value",
                "Percent query coverage per hsp",
                CArgDescriptions::eDouble);
    arg_desc.SetConstraint(kArgQueryCovHspPerc, CArgAllow_Doubles(0.0, 100.0));

    arg_desc.AddOptionalKey(kArgQueryCovHspRes, "int_value",
                "Residues query coverage per hsp",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgQueryCovHspRes, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value",
                "Set maximum number of HSPs per subject sequence to save for each query",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMaxHSPsPerSubject, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddOptionalKey(kArgMaxTargetSequences, "num_sequences",
                "Maximum number of aligned sequences to keep \n"
                "(value of 5 or more is recommanded)\n"
                "Default = '" + NStr::IntToString(BLAST_HITLIST_SIZE) + "'",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMaxTargetSequences, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddFlag(kArgSubjectBestHit, "Turn on best hit per subject sequence", true);

    /// statistical options
    arg_desc.SetCurrentGroup(kGroupStatistics);

    arg_desc.AddOptionalKey(kArgEffSearchSpace, "int_value",
                "Effective length of the search space",
                CArgDescriptions::eInt8);
    arg_desc.SetConstraint(kArgEffSearchSpace, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgDbSize, "num_letters",
                "Effective length of the database",
                CArgDescriptions::eInt8);

    /// extension options
    arg_desc.SetCurrentGroup(kGroupExtension);
    
    arg_desc.AddFlag(kArgSkipOverhang, "Skip overhangs in alignments", true);    

    /// miscellaneous options
    arg_desc.SetCurrentGroup(kGroupMiscellaneous);

    arg_desc.AddDefaultKey(kArgOutput, "output_path",
                "results file",
                CArgDescriptions::eString,
                kDfltOutput);

    arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
                "Number of threads (CPUs) to use in the search",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltNumThreads));
    arg_desc.SetConstraint(kArgNumThreads, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddOptionalKey(kArgGrid, "GRID_options",
                "Index of this node when multiple computation nodes are used in the search\n"
                "(Format: 'node_id num_nodes')\n"
                "Default = '0 1'",
                CArgDescriptions::eString);
}

void CommandLineArguments::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& options)
{
    /// general search options
    if (args.Exist(kTask) && args[kTask].HasValue()) {
        m_Options->task = ProgramNameToEnum(args[kTask].AsString());
    }

    if (args.Exist(kArgEvalue) && args[kArgEvalue].HasValue()) {
        m_Options->expect_value = args[kArgEvalue].AsDouble();
    }

    if (args.Exist(kArgMismatch) && args[kArgMismatch].HasValue()) {
        m_Options->penalty = args[kArgMismatch].AsInteger();
    }

    if (args.Exist(kArgMatch) && args[kArgMatch].HasValue()) {
        m_Options->reward = args[kArgMatch].AsInteger();
    }

    if (args.Exist(kArgGapOpen) && args[kArgGapOpen].HasValue()) {
        m_Options->gap_open = args[kArgGapOpen].AsInteger();
    }

    if (args.Exist(kArgGapExtend) && args[kArgGapExtend].HasValue()) {
        m_Options->gap_extend = args[kArgGapExtend].AsInteger();
    }

    if (args.Exist(kArgSkipOverhang))
        m_Options->skip_overhang = static_cast<bool>(args[kArgSkipOverhang]);

    /// input query options
    if (args.Exist(kArgStrand) && args[kArgStrand].HasValue()) {
        string strand = args[kArgStrand].AsString();
        NStr::ToUpper(strand);
        if (strand == "BOTH") {
            m_Options->strand = F_R;
        } else if (strand == "PLUS") {
            m_Options->strand = FWD;
        } else if (strand == "MINUS") {
            m_Options->strand = REV;
        } else {
            HBN_ERR("Invalid strand value '%s'", strand.c_str());
        }
    }

    /// database options
    if (args.Exist(kArgDbDir) && args[kArgDbDir].HasValue()) {
        if (m_Options->db_dir) free((void*)m_Options->db_dir);
        string db_dir = args[kArgDbDir].AsString();
        m_Options->db_dir = strdup(db_dir.c_str());
    }

    if (args.Exist(kArgKeepDb))
        m_Options->keep_db = static_cast<bool>(args[kArgKeepDb]);

    if (args.Exist(kArgMinQuerySize) && args[kArgMinQuerySize].HasValue()) {
        m_Options->min_query_size = args[kArgMinQuerySize].AsInteger();
    }

    if (args.Exist(kArgMaxQueryVolSeqs) && args[kArgMaxQueryVolSeqs].HasValue()) {
        m_Options->max_query_vol_seqs = args[kArgMaxQueryVolSeqs].AsInteger();
    }

    if (args.Exist(kArgMaxQueryVolRes) && args[kArgMaxQueryVolRes].HasValue()) {
        m_Options->max_query_vol_res = args[kArgMaxQueryVolRes].AsInt8();
    }

    if (args.Exist(kArgMinSubjectSize) && args[kArgMinSubjectSize].HasValue()) {
        m_Options->min_subject_size = args[kArgMinSubjectSize].AsInteger();
    }

    if (args.Exist(kArgMaxSubjectVolSeqs) && args[kArgMaxSubjectVolSeqs].HasValue()) {
        m_Options->max_subject_vol_seqs = args[kArgMaxSubjectVolSeqs].AsInteger();
    }

    if (args.Exist(kArgMaxSubjectVolRes) && args[kArgMaxSubjectVolRes].HasValue()) {
        m_Options->max_subject_vol_res = args[kArgMaxSubjectVolRes].AsInt8();
    }

    /// ddf scoring options
    if (args.Exist(kArgKmerSize) && args[kArgKmerSize].HasValue()) {
        m_Options->kmer_size = args[kArgKmerSize].AsInteger();
    }    

    if (args.Exist(kArgKmerWindow) && args[kArgKmerWindow].HasValue()) {
        m_Options->kmer_window = args[kArgKmerWindow].AsInteger();
    }

    if (args.Exist(kArgBlockSize) && args[kArgBlockSize].HasValue()) {
        m_Options->block_size = args[kArgBlockSize].AsInteger();
    }

    if (args.Exist(kArgDDFScore) && args[kArgDDFScore].HasValue()) {
        m_Options->ddf_score = args[kArgDDFScore].AsInteger();
    }

    if (args.Exist(kArgMaxKmerOcc) && args[kArgMaxKmerOcc].HasValue()) {
        m_Options->max_kmer_occ = args[kArgMaxKmerOcc].AsInteger();
    }

    /// mem chaining scoring options
    if (args.Exist(kArgMemScKmerSize) && args[kArgMemScKmerSize].HasValue()) {
        m_Options->memsc_kmer_size = args[kArgMemScKmerSize].AsInteger();
    }

    if (args.Exist(kArgMemScKmerWindow) && args[kArgMemScKmerWindow].HasValue()) {
        m_Options->memsc_kmer_window = args[kArgMemScKmerWindow].AsInteger();
    }

    if (args.Exist(kArgMemScMemScore) && args[kArgMemScMemScore].HasValue()) {
        m_Options->memsc_score = args[kArgMemScMemScore].AsInteger();
    }

    /// output format
    if (args.Exist(kArgOutputFormat) && args[kArgOutputFormat].HasValue()) {
        string fmtstr = args[kArgOutputFormat].AsString();
        int fmt = NStr::StringToInt(fmtstr);
        switch (fmt) {
            case ePairwise:
                m_Options->outfmt = ePairwise;
                break;
            case eTabular:
                m_Options->outfmt = eTabular;
                break;
            case eTabularWithComments:
                m_Options->outfmt = eTabularWithComments;
                break;
            case eSAM:
                m_Options->outfmt = eSAM;
                break;
            default:
                HBN_ERR("unsupported output format '%d'", fmt);
                break;
        }
    }

    if (args.Exist(kArgSamRgInfo) && args[kArgSamRgInfo].HasValue()) {
        m_Options->rg_info = strdup(args[kArgSamRgInfo].AsString().c_str());
    }

    if (args.Exist(kArgSamRgSample) && args[kArgSamRgSample].HasValue()) {
        m_Options->rg_sample = strdup(args[kArgSamRgSample].AsString().c_str());
    }

    if (args.Exist(kArgSamDumpMd))
        m_Options->dump_md = static_cast<bool>(args[kArgSamDumpMd]);

    /// query filtering options
    if (args.Exist(kArgDustFiltering) && args[kArgDustFiltering].HasValue()) {
        string duststr = args[kArgDustFiltering].AsString();
        CTempString cduststr(duststr);
        CTempString delim(" ");
        vector<string> dust_components;
        NStr::Split(cduststr, delim, dust_components);
        if (dust_components.size() == 1) {
            string s = dust_components[0];
            NStr::ToUpper(s);
            if (s == "NO") {
                m_Options->use_dust_masker = 0;
            } else if (s == "YES") {

            } else {
                HBN_ERR("Invalid value '%s' to argument '%s'", duststr.c_str(), kArgDustFiltering.c_str());
            }
        } else if (dust_components.size() == 3) {
            m_Options->dust_level = NStr::StringToInt(dust_components[0]);
            m_Options->dust_window = NStr::StringToInt(dust_components[1]);
            m_Options->dust_linker = NStr::StringToInt(dust_components[2]);
        } else {
            HBN_ERR("Invalid value '%s' to argument '%s'", duststr.c_str(), kArgDustFiltering.c_str());
        }
    }

    if (args.Exist(kArgUseLCaseMasking))
        m_Options->use_lower_case_masker = static_cast<bool>(args[kArgUseLCaseMasking]);

    /// restrict search or results
    if (args.Exist(kArgPercentIdentity) && args[kArgPercentIdentity].HasValue()) {
        m_Options->perc_identity = args[kArgPercentIdentity].AsDouble();
    }

    if (args.Exist(kArgQueryCovHspPerc) && args[kArgQueryCovHspPerc].HasValue()) {
        m_Options->query_cov_hsp_perc = args[kArgQueryCovHspPerc].AsDouble();
    }

    if (args.Exist(kArgQueryCovHspRes) && args[kArgQueryCovHspRes].HasValue()) {
        m_Options->query_cov_hsp_res = args[kArgQueryCovHspRes].AsInteger();
    }

    if (args.Exist(kArgMaxHSPsPerSubject) && args[kArgMaxHSPsPerSubject].HasValue()) {
        m_Options->max_hsps_per_subject = args[kArgMaxHSPsPerSubject].AsInteger();
    }

    if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences].HasValue()) {
        m_Options->hitlist_size = args[kArgMaxTargetSequences].AsInteger();
    }

    if (args.Exist(kArgSubjectBestHit))
        m_Options->keep_best_hsp_per_subject = static_cast<bool>(args[kArgSubjectBestHit]);

    /// statistical options
    if (args.Exist(kArgEffSearchSpace) && args[kArgEffSearchSpace].HasValue()) {
        m_Options->searchsp_eff = args[kArgEffSearchSpace].AsInt8();
    }

    if (args.Exist(kArgDbSize) && args[kArgDbSize].HasValue()) {
        m_Options->db_length = args[kArgDbSize].AsInt8();
    }

    /// misc options

    if (args.Exist(kArgOutput) && args[kArgOutput].HasValue()) {
        m_Options->output = strdup(args[kArgOutput].AsString().c_str());
    }

    const int kMaxValue = static_cast<int>(hbn_get_cpu_count());
    if (args.Exist(kArgNumThreads) &&
        args[kArgNumThreads].HasValue()) {  // could be cancelled by the exclusion in CRemoteArgs

        // use the minimum of the two: user requested number of threads and
        // number of available CPUs for number of threads
        int num_threads = args[kArgNumThreads].AsInteger();
        if (num_threads > kMaxValue) {
            m_Options->num_threads = kMaxValue;

            string warn_msg = (string)"Number of threads was reduced to " +
                     NStr::IntToString((unsigned int)m_Options->num_threads) +
                     " to match the number of available CPUs";
            HBN_WARN("%s", warn_msg.c_str());
        }
        else {
            m_Options->num_threads = num_threads;
        }
    }

    if (args.Exist(kArgGrid) && args[kArgGrid].HasValue()) {
        string gridstr = args[kArgGrid].AsString();
        CTempString kDelim(" ");
        vector<string> components;
        NStr::Split(gridstr, " ", components);
        if (components.size() != 2) 
            HBN_ERR("Invalid value '%s' to argument '%s'", gridstr.c_str(), kArgGrid.c_str());
        m_Options->node_id = NStr::StringToInt(components[0]);
        m_Options->num_nodes = NStr::StringToInt(components[1]);

        if (m_Options->node_id < 0) HBN_ERR("node id must be >=0: %s", gridstr.c_str());
        if (m_Options->num_nodes <= 0) HBN_ERR("number of nodes must be >0: %s", gridstr.c_str());
        if (m_Options->node_id >= m_Options->num_nodes) 
            HBN_ERR("node index (%d) must be smaller than number of nodes (%d)", 
                m_Options->node_id, m_Options->num_nodes);
    }
}

void Init_HbnProgramOptions(HbnProgramOptions* opts)
{
    /// general search options
    opts->task = kDfltTask;
    opts->expect_value = kDfltEvalue;
    opts->penalty = kDfltMismatch;
    opts->reward = kDfltMatch;
    opts->gap_open = kDfltGapOpen;
    opts->gap_extend = kDfltGapExtend;
    opts->skip_overhang = kDfltSkipOverhang;

    /// input query options
    opts->strand = F_R;

    /// database options
    opts->db_dir = strdup(kDfltDbDir.c_str());
    opts->keep_db = kDfltKeepDb;
    opts->min_query_size = kDfltMinQuerySize;
    opts->max_query_vol_seqs = kDfltMaxQueryVolSeqs;
    opts->max_query_vol_res = kDfltMaxQueryVolRes;
    opts->min_subject_size = kDfltMinSubjectSize;
    opts->max_subject_vol_seqs = kDfltMaxSubjectVolSeqs;
    opts->max_subject_vol_res = kDfltMaxSubjectVolRes;

    /// ddf scoring options
    opts->kmer_size = kDfltKmerSize;
    opts->kmer_window = kDfltKmerWindow;
    opts->max_kmer_occ = kDfltMaxKmerOcc;
    opts->block_size = kDfltBlockSize;
    opts->ddf_score = kDfltDDFScore;

    /// mem chaining scoring options
    opts->memsc_kmer_size = kDfltMemScKmerSize;
    opts->memsc_kmer_window = kDfltMemScKmerWindow;
    opts->memsc_score = kDfltMemScMemScore;

    /// formatting options
    opts->outfmt = kDfltOutfmt;
    opts->rg_info = kDfltSamRgInfo;
    opts->rg_sample = kDfltSamRgSample;
    opts->dump_md = kDfltSamDumpMd;

    /// query filtering options
    opts->use_dust_masker = 1;
    opts->dust_level = kDustLevel;
    opts->dust_linker = kDustLinker;
    opts->dust_window = kDustWindow;
    opts->use_lower_case_masker = kDfltArgUseLCaseMasking;

    /// restrict search or results
    opts->perc_identity = kDfltPercentIdentity;
    opts->query_cov_hsp_perc = kDfltQueryCovHspPerc;
    opts->query_cov_hsp_res = kDfltQueryCovHspRes;
    opts->max_hsps_per_subject = kDfltMaxHSPsPerSubject;
    opts->hitlist_size = kDfltMaxTargetSequences;
    opts->keep_best_hsp_per_subject = kDfltSubjectBestHit;

    /// statistical options
    opts->searchsp_eff = 0;
    opts->db_length = 0;
    opts->dbseq_num = 0;

    /// misc options
    opts->num_threads = kDfltNumThreads;
    opts->node_id = kDfltNodeId;
    opts->num_nodes = kDfltNumNodes;

    opts->query = NULL;
    opts->subject = NULL;
    opts->output = kDfltOutput.c_str();
}

static void
s_PreCheckCmdLineArgs(int argc, char* argv[], CArgDescriptions* arg_desc)
{
    string kProgram = FindProgramDisplayName(argv[0]);
    string kHbnUsage = kProgram + " [OPTIONS] query subject";

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') continue;
        if (NStr::CompareCase(argv[i] + 1, kArgHelp) == 0) {
            string usage_info;
            //arg_desc->PrintUsage(usage_info);
            arg_desc->HbnPrintUsage(kHbnUsage, usage_info);
            cout << usage_info << endl;
            exit(0);
        } else if (NStr::CompareCase(argv[i] + 1, kArgFullHelp) == 0) {
            string usage_info;
            //arg_desc->PrintUsage(usage_info, true);
            arg_desc->HbnPrintUsage(kHbnUsage, usage_info, true);
            cout << usage_info << endl;
            exit(0);
        } else if (NStr::CompareCase(argv[i] + 1, kArgVersion) == 0) {
            string appname = FindProgramDisplayName(argv[0]);
            string version = PrintProgramVersion(appname);
            cout << version << endl;
            exit(0);
        }
    }
}

extern "C"
void ParseHbnProgramCmdLineArguments(int argc, char* argv[], HbnProgramOptions* opts)
{
    Init_HbnProgramOptions(opts);
    TBlastCmdLineArgs arg_list;

    /// setup description
    CRef<IBlastCmdLineArgs> arg;
    string kProgram = FindProgramDisplayName(argv[0]);
    string kProgramDescription("Nucleotide-Nucleotide sequence alignment toolkit");
    arg.reset(new CProgramDescriptionArgs(kProgram, kProgramDescription));
    arg_list.push_back(arg);

    arg.reset(new CommandLineArguments(opts));
    arg_list.push_back(arg);

    unique_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);
    NON_CONST_ITERATE(TBlastCmdLineArgs, arg_iter, arg_list) {
        (*arg_iter)->SetArgumentDescriptions(*arg_desc);
    }

    /// examine trivial arguments (-help, -h, -version)
    s_PreCheckCmdLineArgs(argc, argv, arg_desc.get());

    /// process nontrivial arguments
    unique_ptr<CArgs> cmd_args(new CArgs);
    int argv_idx = 1;
    auto& supported_args = (*arg_desc).GetArgs();
    while (argv_idx < argc) {
        if (argv[argv_idx][0] != '-') break;
        string argname = argv[argv_idx] + 1;
        //cout << "process " << argname << endl;
        //if (!(*arg_desc).Exist(argname)) HBN_ERR("unrecognised argument '%s'", argname.c_str());
        bool negative = false;
        auto it = (*arg_desc).Find(argname, &negative);
        if (it == supported_args.end()) HBN_ERR("unrecognised argument '%s'", argname.c_str());
        CArgDesc& arg = **it;
        CArgValue* av = nullptr;
        if (ArgDescIsFlag(arg)) {
            //av = arg.ProcessDefault();
            av = arg.ProcessArgument(kEmptyStr);
            ++argv_idx;
        } else {
            if (argv_idx + 1 >= argc) HBN_ERR("Mandatory value to argument '%s' is missing", argname.c_str());
            av = arg.ProcessArgument(argv[argv_idx + 1]);
            argv_idx += 2;
        }
        cmd_args->Add(av, true, true);
    }

    CBlastOptions cblastopts;
    NON_CONST_ITERATE(TBlastCmdLineArgs, arg_iter, arg_list) {
        (*arg_iter)->ExtractAlgorithmOptions(*cmd_args, cblastopts);
    }

    /// query, subject
    if (argc - argv_idx  < 2) {
        HBN_ERR("The query and subject must be specified");
    } else if (argc - argv_idx > 2) {
        string err = "Too many query and subject values: '";
        for (int i = argv_idx; i < argc; ++i) {
            err += argv[i];
            if (i != argc - 1) err += ' ';
        }
        err += "'";
        HBN_ERR("%s", err.c_str());
    } else {
        opts->query = argv[argv_idx];
        opts->subject = argv[argv_idx + 1];
    }
}

#define os_one_option_value(name, value) os << '-' << name << ' ' << value << ' '
#define os_one_flag_option(name) os << '-' << name << ' ';

extern "C"
char* HbnProgramOptions2String(const HbnProgramOptions* opts)
{
    ostringstream os;

    /// general search options
    os_one_option_value(kTask, EProgramToTaskName(opts->task));
    os_one_option_value(kArgEvalue, opts->expect_value);
    os_one_option_value(kArgMismatch, opts->penalty);
    os_one_option_value(kArgMatch, opts->reward);
    os_one_option_value(kArgGapOpen, opts->gap_open);
    os_one_option_value(kArgGapExtend, opts->gap_extend);
    if (opts->skip_overhang) os_one_flag_option(kArgSkipOverhang);

    /// input query options
    const char* strand_name[] = { "plus", "minus", "both" };
    os_one_option_value(kArgStrand, strand_name[opts->strand]);

    /// database options
    os_one_option_value(kArgDbDir, opts->db_dir);
    if (opts->keep_db) os_one_flag_option(kArgKeepDb);
    if (opts->min_query_size) os_one_option_value(kArgMinQuerySize, opts->min_query_size);
    if (opts->max_query_vol_seqs != kDfltMaxQueryVolSeqs) os_one_option_value(kArgMaxQueryVolSeqs, opts->max_query_vol_seqs);
    string size_str = NStr::UInt8ToString_DataSize(opts->max_query_vol_res);
    os_one_option_value(kArgMaxQueryVolRes, size_str);
    if (opts->min_subject_size) os_one_option_value(kArgMinSubjectSize, opts->min_subject_size);
    if (opts->max_subject_vol_seqs != kDfltMaxSubjectVolSeqs) os_one_option_value(kArgMaxSubjectVolSeqs, opts->max_subject_vol_seqs);
    size_str = NStr::UInt8ToString_DataSize(opts->max_subject_vol_res);
    os_one_option_value(kArgMaxSubjectVolRes, size_str);

    /// ddf scoring
    os_one_option_value(kArgKmerSize, opts->kmer_size);
    os_one_option_value(kArgKmerWindow, opts->kmer_window);
    os_one_option_value(kArgMaxKmerOcc, opts->max_kmer_occ);
    os_one_option_value(kArgBlockSize, opts->block_size);
    os_one_option_value(kArgDDFScore, opts->ddf_score);

    /// mem chaining scoring options
    os_one_option_value(kArgMemScKmerSize, opts->memsc_kmer_size);
    os_one_option_value(kArgMemScKmerWindow, opts->memsc_kmer_window);
    os_one_option_value(kArgMemScMemScore, opts->memsc_score);

    /// output format
    os_one_option_value(kArgOutputFormat, opts->outfmt);

    /// query filtering options
    if (opts->use_dust_masker) {
        os << '-' << kArgDustFiltering << ' '
           << "'" << opts->dust_level 
           << " " << opts->dust_window
           << " " << opts->dust_linker
           << "' ";
    }
    if (opts->use_lower_case_masker) os_one_flag_option(kArgUseLCaseMasking);

    /// restrict search of results
    os_one_option_value(kArgPercentIdentity, opts->perc_identity);
    os_one_option_value(kArgQueryCovHspPerc, opts->query_cov_hsp_perc);
    os_one_option_value(kArgQueryCovHspRes, opts->query_cov_hsp_res);
    os_one_option_value(kArgMaxHSPsPerSubject, opts->max_hsps_per_subject);
    os_one_option_value(kArgMaxTargetSequences, opts->hitlist_size);
    if (opts->keep_best_hsp_per_subject) os_one_flag_option(kArgSubjectBestHit);

    /// misc options
    os_one_option_value(kArgNumThreads, opts->num_threads);
    os << '-' << kArgGrid << ' ' << opts->node_id << ' ' << opts->num_nodes << ' ';

    size_str = os.str();
    return strdup(size_str.c_str());
}

END_NCBI_SCOPE