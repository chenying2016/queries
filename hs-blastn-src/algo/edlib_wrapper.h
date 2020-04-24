#ifndef __EDLIB_WRAPPER_H
#define __EDLIB_WRAPPER_H

#include "../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    kstring_t       query;
    kstring_t       target;
    kstring_t       qaln;
    kstring_t       taln;
    int             tolerance;
    int             do_traceback;
    int             dist;
    int             score;
} EdlibAlignData;

EdlibAlignData*
EdlibAlignDataNew();

EdlibAlignData*
EdlibAlignDataFree(EdlibAlignData* data);

int
edlib_nw(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    kstring_t* qaln,
    kstring_t* taln);

int
edlib_shw(EdlibAlignData* data,
    const u8* query,
    int query_size,
    const u8* target,
    int target_size,
    int* qend,
    int* tend,
    kstring_t* qaln,
    kstring_t* taln);

void
edlib_extend(EdlibAlignData* data,
    const u8* query,
    const int query_size,
    const u8* target,
    const int target_size,
    const int block_size,
    const BOOL right_extend,
    vec_u8* qfrag,
    vec_u8* tfrag,
    int* qend,
    int* tend,
    kstring_t* qaln,
    kstring_t* taln);

#ifdef __cplusplus
}
#endif

#endif // __EDLIB_WRAPPER_H