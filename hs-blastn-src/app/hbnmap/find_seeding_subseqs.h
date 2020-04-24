#ifndef __FIND_SEEDING_SUBSEQS_H
#define __FIND_SEEDING_SUBSEQS_H

#include "../../corelib/hbn_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

int find_seeding_subseqs(const unsigned char* seq, 
        const unsigned int seq_size,
        const int min_size,
        vec_int_pair* good_regions);

#ifdef __cplusplus
}
#endif

#endif // __FIND_SEEDING_SUBSEQS_H