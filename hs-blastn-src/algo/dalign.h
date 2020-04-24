#ifndef __DALIGN_H
#define __DALIGN_H

#include "align.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    Align_Spec*     align_spec;
    Work_Data*      work;
    float           basis_freq[4];
    Path            path;
    Alignment       align;
    kstring_t       aseq;
    kstring_t       bseq;
    double          error;
    double          ident_perc;
} DalignData;

#define ocda_ident_perc(ocda)     ((ocda).ident_perc)
#define ocda_query_start(ocda)    ((ocda).path.abpos)
#define ocda_query_end(ocda)      ((ocda).path.aepos)
#define ocda_target_start(ocda)   ((ocda).path.bbpos)
#define ocda_target_end(ocda)     ((ocda).path.bepos)
#define ocda_distance(ocda)		  ((ocda).path.diffs)

DalignData*
DalignDataFree(DalignData* data);

DalignData*
DalignDataNew(double error);

int
dalign_align(DalignData* ocda,
	const u8* query,
    int qoff,
    int qsize,
	const u8* target,
    int toff,
    int tsize,
	const int min_align_size,
	const double min_ident_perc,
	int* qbeg,
	int* qend,
	int* tbeg,
	int* tend,
	double* ident_perc,
	kstring_t* qaln,
	kstring_t* taln);

#ifdef __cplusplus
}
#endif

#endif // __DALIGN_H