#ifndef __HBN_OPTIONS_HANDLE_H
#define __HBN_OPTIONS_HANDLE_H

#include "hbn_options.h"
#include "../../ncbi_blast/setup/blast_options.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    EProgram                        m_Program;
    QuerySetUpOptions*              m_QueryOpts;
    LookupTableOptions*             m_LutOpts;
    BlastInitialWordOptions*        m_InitWordOpts;
    BlastExtensionOptions*          m_ExtnOpts;
    BlastHitSavingOptions*          m_HitSaveOpts;
    PSIBlastOptions*                m_PSIBlastOpts;
    BlastDatabaseOptions*           m_DbOpts;
    BlastScoringOptions*            m_ScoringOpts;
    BlastEffectiveLengthsOptions*   m_EffLenOpts;
} HbnOptionsHandle;

HbnOptionsHandle*
HbnOptionsHandleNew(EProgram program);

HbnOptionsHandle*
HbnOptionsHandleFree(HbnOptionsHandle* opts);

void
HbnOptionsHandle_Update(const HbnProgramOptions* opts, HbnOptionsHandle* handle);

#ifdef __cplusplus
}
#endif

#endif // __HBN_OPTIONS_HANDLE_H