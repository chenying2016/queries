#ifndef __PRIMER_MAP_ONE_VOLUME_H
#define __PRIMER_MAP_ONE_VOLUME_H

#include "../hbnmap/hbn_options.h"
#include "../hbnmap/hbn_options_handle.h"
#include "../../corelib/seqdb.h"

#ifdef __cplusplus
extern "C" {
#endif

void
primer_map_one_volume(CSeqDB* qvol, 
    CSeqDB* svol, 
    const HbnProgramOptions* opts,
    const HbnOptionsHandle* opts_handle,
    FILE* out);

#ifdef __cplusplus
}
#endif

#endif // __PRIMER_MAP_ONE_VOLUME_H