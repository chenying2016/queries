#ifndef __HBN_TRACEBACK_AUX_H
#define __HBN_TRACEBACK_AUX_H

#include "../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define dump_align_string(qaln, saln, aln_size, stream) do { \
	hbn_fwrite(qaln, 1, aln_size, stream); \
	fprintf(stream, "\n"); \
	for (int _i = 0; _i < aln_size; ++_i) { \
		if (qaln[_i] == saln[_i]) fprintf(stream, "|"); \
		else fprintf(stream, "*"); \
	} \
	fprintf(stream, "\n"); \
	hbn_fwrite(saln, 1, aln_size, stream); \
	fprintf(stream, "\n"); \
} while (0)

double
calc_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size,
                int* dist,
				int* score);

void
validate_aligned_string(const char* source_file,
						const char* source_func,
						const int source_line,
						int qid,
						const u8* query,
						const int qoff,
						const int qend,
						const char* query_mapped_string,
						int tid,
						const u8* target,
						const int toff,
						const int tend,
						const char* target_mapped_string,
						const size_t align_size,
					    const BOOL right_extend);

#ifdef __cplusplus
}
#endif

#endif // __HBN_TRACEBACK_AUX_H