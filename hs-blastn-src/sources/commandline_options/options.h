#ifndef BLAST_OPTIONS_H
#define	BLAST_OPTIONS_H

#include "def.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>

/* filtering options */

struct SeedingOptions
{
    int seed_size;
};

struct PrintVersionOption
{
    Boolean print_version;
};

struct HelpOptions
{
    Boolean simple_help;
    Boolean full_help;
};

struct InputOptions
{
    const char* query;
	const char* query_list;
    const char* db;
};

/// Defines the output formats supported by our command line formatter
enum EOutputFormat {
	/// Standard pairwise alignments
	ePairwise = 0,                  
	///< Query anchored showing identities
	eQueryAnchoredIdentities,
	///< Query anchored no identities
	eQueryAnchoredNoIdentities,
	///< Flat query anchored showing identities
	eFlatQueryAnchoredIdentities,
	///< Flat query anchored no identities
	eFlatQueryAnchoredNoIdentities,
	/// XML output
	eXml,
	/// Tabular output
	eTabular,
	/// Tabular output with comments
	eTabularWithComments,
	/// ASN.1 text output
	eAsnText,
	/// ASN.1 binary output
	eAsnBinary,
	/// Comma-separated values
	eCommaSeparatedValues,
	/// BLAST archive format
	eArchiveFormat,
	/// JSON seq-align
	eJsonSeqalign,
	/// JSON XInclude
	eJson,
	/// XML2 XInclude
	eXml2,
	/// JSON2 single file
	eJson_S,
	/// XML2 single file
	eXml2_S,
	/// SAM format
	eSAM,
	/// Sentinel value for error checking
	eEndValue
};

struct OutputOptions
{
    const char* output_file_name;
    FILE* out;
    EOutputFormat outfmt;
    int num_descriptions;
    int num_alignments;
};

enum ESequenceStrand
{
    eStrandPlus = 0,
    eStrandMinus = 1,
    eStrandBoth = 2
};

struct RunningOptions
{
    int num_threads;
    ESequenceStrand query_strand;
};

struct SDustOptions {
    int level;
    int window;
    int linker;

    SDustOptions(int level_ = kDustLevel, int window_ = kDustWindow, int linker_ = kDustLinker)
    : level(level_), window(window_), linker(linker_) {
    }

    static const int kDustLevel = 20;
    static const int kDustWindow = 64;
    static const int kDustLinker = 1;
};

#define kDefaultRepeatFilterDb "repeat/repeat_9606"

struct SRepeatFilterOptions {
    char* database;

    SRepeatFilterOptions() {
        database = strdup(kDefaultRepeatFilterDb);
    }

    ~SRepeatFilterOptions() {
        sfree(database);
    }

    Int2 ResetDB(const char* db) {
        database = strdup(db);
        return 0;
    }
};

struct SWindowMaskerOptions {
    int taxid; // select masking database for this taxid
    const char* database; // use winmasker database at this location

    SWindowMaskerOptions(int tid = 0, const char* db = NULL)
    : taxid(tid), database(db) {
    }

    ~SWindowMaskerOptions() {
        if (database) database = NULL; //sfree(database);
    }

    Int2 ResetDB(const char* db) {
        sfree(database);
        if (db)
            database = strdup(db);
        return 0;
    }
};

struct SBlastFilterOptions {
    Boolean mask_at_seeding;
    SDustOptions* dustOptions;
    SRepeatFilterOptions* repeatFilterOptions;
    SWindowMaskerOptions* windowMaskerOptions;

    SBlastFilterOptions() : mask_at_seeding(TRUE) {
        dustOptions = new SDustOptions();
        repeatFilterOptions = new SRepeatFilterOptions();
        windowMaskerOptions = new SWindowMaskerOptions();
    }

    ~SBlastFilterOptions() {
        delete dustOptions;
        delete repeatFilterOptions;
        delete windowMaskerOptions;
    }
};

struct BlastScoringOptions {
    Int2 reward; /**< Reward for a match */
    Int2 penalty; /**< Penalty for a mismatch */
    Boolean gapped_calculation; /**< gap-free search if FALSE */
    Int4 gap_open; /**< Extra penalty for starting a gap */
    Int4 gap_extend; /**< Penalty for each gap residue */
};

/** Options needed for initial word finding and processing */
struct BlastInitialWordOptions {
    double gap_trigger; /**< Score in bits for starting gapped extension */
    Int4 window_size; /**< Maximal allowed distance between 2 hits in case 2
                        hits are required to trigger the extension */
    Int4 scan_range; /**< Maximal number of gaps allowed between 2 hits */
    double x_dropoff; /**< X-dropoff value (in bits) for the ungapped
                         extension */
};

/** The algorithm to be used for preliminary
 *  gapped extensions
 */
enum EBlastPrelimGapExt {
    eDynProgScoreOnly, /**< standard affine gapping */
    eGreedyScoreOnly, /**< Greedy extension (megaBlast) */
    eSmithWatermanScoreOnly /**< Score-only smith-waterman */
};

/** The algorithm to be used for final gapped
 *  extensions with traceback
 */
enum EBlastTbackExt {
    eDynProgTbck, /**< standard affine gapping */
    eGreedyTbck, /**< Greedy extension (megaBlast) */
    eSmithWatermanTbck, /**< Smith-waterman finds optimal scores, then
                                ALIGN_EX to find alignment. */
    eSmithWatermanTbckFull /**< Smith-waterman to find all alignments */
};

/** Options used for gapped extension
 *  These include:
 *  a. Penalties for various types of gapping;
 *  b. Drop-off values for the extension algorithms tree exploration;
 *  c. Parameters identifying what kind of extension algorithm(s) should
 *     be used.
 */
struct BlastExtensionOptions {
    double gap_x_dropoff; /**< X-dropoff value for gapped extension (in bits) */
    double gap_x_dropoff_final; /**< X-dropoff value for the final gapped
                                  extension (in bits) */
    EBlastPrelimGapExt ePrelimGapExt; /**< type of preliminary gapped extension (normally) for calculating
                              score. */
    EBlastTbackExt eTbackExt; /**< type of traceback extension. */
    Int4 compositionBasedStats; /**< mode of compositional adjustment to use;
                                   if zero then compositional adjustment is
                                   not used */
    Int4 unifiedP; /**< Indicates unified P values to be used in blastp or tblastn */
};

/** Options used when evaluating and saving hits
 *  These include:
 *  a. Restrictions on the number of hits to be saved;
 *  b. Restrictions on the quality and positions of hits to be saved;
 *  c. Parameters used to evaluate the quality of hits.
 */
struct BlastHitSavingOptions {
    double expect_value; /**< The expect value cut-off threshold for an HSP, or
                            a combined hit if sum statistics is used */
    Int4 cutoff_score; /**< The (raw) score cut-off threshold */
    double percent_identity; /**< The percent identity cut-off threshold */

    Int4 hitlist_size; /**< Maximal number of database sequences to return
                        results for */
    Int4 hsp_num_max; /**< Maximal number of HSPs to save for one database
                        sequence */
    Int4 total_hsp_limit; /**< Maximal total number of HSPs to keep */
    Int4 culling_limit; /**< If the query range of an HSP is contained in
                            at least this many higher-scoring HSPs, throw
                            away the HSP as redundant (turned off if zero) */
    Int4 mask_level; /**< Only keep the highest scoring HSP when more than
                          one HSP overlaps the same region of the query by
                          more than or equal to mask_level %. -RMH- */

    /********************************************************************/
    /* Merge all these in a structure for clarity? */
    Boolean do_sum_stats; /**< Force sum statistics to be used to combine HSPs,
                          TRUE by default for all ungapped searches and translated
                          gapped searches (except RPS-BLAST) */
    Int4 longest_intron; /**< The longest distance between HSPs allowed for
                           combining via sum statistics with uneven gaps */
    /********************************************************************/

    Int4 min_hit_length; /**< optional minimum alignment length; alignments
                                not at least this long are discarded */
    Int4 min_diag_separation; /**< How many diagonals separate a hit from a substantial alignment
                                  before it's not blocked out. Must be > 0 to be used. */
    //EBlastProgramType program_number; /**< indicates blastn, blastp, etc. */

    /** Contains options to configure the HSP filtering/writering structures
     * If not set, the default HSP filtering strategy is used.
     */
    //BlastHSPFilteringOptions* hsp_filt_opt;

    /** Low-score option.  Do not pass ungapped alignments on for later processing if
     * the hitlist is already full of other alignments unless the ungapped aligment
     * is above the fraction X of the least significant database match.
     * zero should turn this off.
     */
    double low_score_perc;

};

/** Options for setting up effective lengths and search spaces.
 * The values are those the user has specified to override the real sizes.
 */
struct BlastEffectiveLengthsOptions {
    Int8 db_length; /**< Database length to be used for statistical
                         calculations */
    Int4 dbseq_num; /**< Number of database sequences to be used for
                           statistical calculations */
    Int4 num_searchspaces; /**< Number of elements in searchsp_eff, this must be
                            equal to the number of contexts in the search */
    Int8 *searchsp_eff; /**< Search space to be used for statistical
                           calculations (one such per query context) */
};

PrintVersionOption*
PrintVersionOptionNew();

HelpOptions*
HelpOptionsNew();

SeedingOptions*
SeedingOptionsNew();

InputOptions*
InputOptionsNew();

OutputOptions*
OutputOptionsNew();

RunningOptions*
RunningOptionsNew();

BlastInitialWordOptions*
BlastInitialWordOptionsNew();

BlastExtensionOptions*
BlastExtensionOptionsNew();

BlastHitSavingOptions*
BlastHitSavingOptionsNew();

BlastEffectiveLengthsOptions*
BlastEffectiveLengthsOptionsNew();

BlastScoringOptions*
BlastScoringOptionsNew();

void PrintScoringOptions(BlastScoringOptions* options);

void
PrintInitialWordOptions(BlastInitialWordOptions* options);

void
PrintExtensionOptions(BlastExtensionOptions* options);

void
PrintHitSavingOptions(BlastHitSavingOptions* options);

void
PrintEffectiveLengthOptions(BlastEffectiveLengthsOptions* options);

Int2 VersionOptionsValidate(const PrintVersionOption* version_options);

Int2 HelpOptionsValidate(const HelpOptions* help_options);

Int2 BLAST_ValidateOptions(const BlastExtensionOptions* ext_options,
                           const BlastScoringOptions* score_options,
                           const BlastInitialWordOptions* word_options,
                           const BlastHitSavingOptions* hit_options,
                           const InputOptions* ioptions,
                           const OutputOptions* ooptions,
                           const SeedingOptions* seed_options,
                           const RunningOptions* running_options,
                           const SBlastFilterOptions* filtering_options);

#endif	/* BLAST_OPTIONS_H */

