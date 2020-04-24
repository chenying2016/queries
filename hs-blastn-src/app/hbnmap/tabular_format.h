#ifndef __TABULAR_FORMAT_HPP
#define __TABULAR_FORMAT_HPP

#include "../../corelib/seqdb.h"
#include "../../ncbi_blast/setup/blast_hits.h"
#include "hbn_options.h"

#ifdef __cplusplus
extern "C" {
#endif

void
print_tabular_reports(HbnHSPResults* results, const char* db_name, const CSeqDB* db, const CSeqDB* queries, const EOutputFormat outfmt);

#ifdef __cplusplus
}
#endif

#endif // __TABULAR_FORMAT_HPP