#ifndef __MAP_ONE_VOLUME_H
#define __MAP_ONE_VOLUME_H

#include "../../ncbi_blast/setup/blast_hits.h"
#include "hbn_task_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

void
hbn_align_one_volume(hbn_task_struct* ht_struct);

#ifdef __cplusplus
}
#endif

#endif // __MAP_ONE_VOLUME_H