/* 
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 *
 */

/** @file blast_traceback.h
 * Functions to do gapped alignment with traceback
 */

#ifndef TRACEBACK_H
#define	TRACEBACK_H

#include "sequence.h"
#include "parameters.h"
#include "dbinfo.h"
#include "gapalign.h"
#include "interval_tree.h"
#include "memallocator.h"
#include "utility.h"

Int2
Blast_TraceBackFromOneContextHSPList(HSP** hsp_list,
									 Int4& num_hsps,
									 SmallObjAllocator& soa,
									 QueryInfo* query_info,
									 BlastScoreBlk* sbp,
									 DbInfo* dbinfo,
									 IntervalTree& itree,
									 GreedyAligner* gapped_aligner,
									 const BlastScoringParameters* score_params,
									 const BlastExtensionParameters* ext_params,
									 const BlastHitSavingParameters* hit_params);

void
AdjustSubjectRange(Int8* subject_offset_ptr, 
				   Int8* subject_length_ptr,
				   Int4 query_offset, 
				   Int4 query_length, 
				   Int8* start_shift);

#endif	/* TRACEBACK_H */

