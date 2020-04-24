#include "hbn_traceback_aux.h"

double
calc_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size,
                int* dist,
				int* score)
{
	if (align_size == 0) return 0.0;
	
	int n = 0;
	for (int i = 0; i < align_size; ++i) {
		if (query_mapped_string[i] == target_mapped_string[i]) ++n;
	}
    if (dist) *dist = align_size - n;
	if (score) {
		int a = abs(MATCH_REWARD);
		int b = -abs(MISMATCH_PENALTY);
		int go = -abs(GAP_OPEN);
		int ge = -abs(GAP_EXTEND);
		int i = 0;
		int s = 0;
		while (i < align_size) {
			int qc = query_mapped_string[i];
			int tc = target_mapped_string[i];
			if (qc != GAP_CHAR && tc != GAP_CHAR) {
				s += (qc == tc) ? a : b;
				++i;
				continue;
			}
			if (qc == GAP_CHAR) {
				int l = 0;
				int j = i;
				for (; j < align_size; ++j) {
					qc = query_mapped_string[j];
					tc = target_mapped_string[j];
					if (qc == GAP_CHAR && tc == GAP_CHAR) continue;
					if (qc != GAP_CHAR) break;
					++l;
				}
				s += (go + ge * l);
				i = j;
				continue;
			}
			hbn_assert(tc == GAP_CHAR);
			{
				int l = 0;
				int j = i;
				for (; j < align_size; ++j) {
					qc = query_mapped_string[j];
					tc = target_mapped_string[j];
					if (qc == GAP_CHAR && tc == GAP_CHAR) continue;
					if (tc != GAP_CHAR) break;
					++l;
				}
				s += (go + ge * l);
				i = j;
				continue;
			}
		}
		*score = s;
	}
	return 100.0 * n / align_size;
}

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
					    const BOOL right_extend)
{
	//return;
	int x = qoff, y = toff;
	for (size_t i = 0; i != align_size; ++i) {
		const char qc = query_mapped_string[i];
		if (qc != GAP_CHAR) {
            int c = right_extend ? query[x] : query[-x];
			const char qc1 = DECODE_RESIDUE(c);
            if (qc != qc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, qc = %c, qc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d, align_size = %lu\n",
					  source_file,
					  source_func,
					  source_line,
					  qid,
					  tid,
					  right_extend,
					  i,
					  x,
					  y,
					  qc,
					  qc1,
					  qoff,
					  qend,
					  toff,
					  tend,
					  align_size);
                abort();
            }		  
			++x;
		}
		const char tc = target_mapped_string[i];
		if (tc != GAP_CHAR) {
            int c = right_extend ? target[y] : target[-y];
			const char tc1 = DECODE_RESIDUE(c);
            if (tc != tc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, tc = %c, tc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d\n",
						  source_func,
						  source_func,
						  source_line,
						  qid,
						  tid,
						  right_extend,
						  i,
						  x,
						  y,
						  tc,
						  tc1,
						  qoff,
						  qend,
						  toff,
						  tend);
                abort();
            }
			++y;
		}
	}
}