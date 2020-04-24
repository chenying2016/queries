#ifndef __HBN_JOB_CONTROL_H
#define __HBN_JOB_CONTROL_H

#include "../../corelib/hbn_aux.h"
#include "../../corelib/seqdb.h"
#include "hbn_options.h"

#ifdef __cplusplus
extern "C" {
#endif

const char*
make_qi_vs_sj_results_path(const char* wrk_dir, const char* stage, const int qi, const int sj, char path[]);

FILE*
open_qi_vs_sj_results_file(const char* wrk_dir, const char* stage, const int qi, const int sj, const char* mode);

BOOL 
qi_vs_sj_is_mapped(const char* wrk_dir, const char* stage, const int qi, const int sj);

void
qi_vs_sj_make_mapped(const char* wrk_dir, const char* stage, const int qi, const int sj);

BOOL
all_vs_sj_is_mapped(const char* wrk_dir, 
    const char* stage, 
    const int qi_start,
    const int num_query_vols,
    const int sj,
    const int node_id,
    const int num_nodes);

void
merge_qi_vs_sj_results(const char* wrk_dir, const char* stage, const int qi, const int sj, const CSeqDB* db, const HbnProgramOptions* opts, FILE* out);

void
merge_all_vs_sj_results(const char* wrk_dir, 
    const char* stage, 
    const int qi_start,
    const int num_query_vols,
    const int sj,
    const int node_id,
    const int num_nodes,
    const HbnProgramOptions* opts,
    FILE* out);

#ifdef __cplusplus
}
#endif

#endif // __HBN_JOB_CONTROL_H