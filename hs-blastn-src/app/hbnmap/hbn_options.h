#ifndef __HBN_OPTIONS_H
#define __HBN_OPTIONS_H

#include "../../ncbi_blast/setup/blast_options.h"
#include "../../ncbi_blast/setup/blast_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define HBN_QUERY_CHUNK_SIZE 20

/// Defines the output formats supported by our command line formatter
typedef enum {
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

    eTaxFormat,
        
    ///igblast AIRR rearrangement, 19
    eAirrRearrangement,
    /// Sentinel value for error checking
    eEndValue     
} EOutputFormat;

typedef struct {
    /// general search options
    EProgram            task;
    double              expect_value;
    int                 penalty;
    int                 reward;
    int                 gap_open;
    int                 gap_extend;
    int                 skip_overhang;

    /// input query options
    int                 strand;

    /// database options
    const char*         db_dir;
    int                 keep_db;
    int                 min_query_size;
    int                 max_query_vol_seqs;
    size_t              max_query_vol_res;
    int                 min_subject_size;
    int                 max_subject_vol_seqs;
    size_t              max_subject_vol_res;

    /// ddf scoring options
    int                 kmer_size;
    int                 kmer_window;
    int                 max_kmer_occ;
    int                 block_size;
    int                 ddf_score;

    /// mem scoring options
    int                 memsc_kmer_size;
    int                 memsc_kmer_window;
    int                 memsc_mem_size;
    int                 memsc_score;

    /// formatting options
    EOutputFormat       outfmt;
    int                 dump_md;
    const char*         rg_info;
    const char*         rg_sample;

    /// query filtering options
    int                 use_dust_masker;
    int                 dust_level;
    int                 dust_window;
    int                 dust_linker;
    int                 use_lower_case_masker;

    /// restrict search or results
    double              perc_identity;
    double              query_cov_hsp_perc;
    int                 query_cov_hsp_res;
    int                 max_hsps_per_subject;
    int                 hitlist_size;
    int                 keep_best_hsp_per_subject;

    /// statistical options
    i64                 searchsp_eff;
    i64                 db_length;
    int                 dbseq_num;

    /// misc options
    int                 num_threads;
    int                 node_id;
    int                 num_nodes;

    const char*         query;
    const char*         subject;
    const char*         output;
} HbnProgramOptions;

#ifdef __cplusplus
}
#endif

#endif // __HBN_OPTIONS_H