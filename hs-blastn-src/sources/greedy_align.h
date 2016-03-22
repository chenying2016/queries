/* ===========================================================================
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
* Author: Ilya Dondoshansky
*
*/

/** @file greedy_align.h
* Prototypes and structures for greedy gapped alignment
*/

#ifndef GREEDY_ALIGN_H
#define	GREEDY_ALIGN_H

#include "def.h"
#include "stat.h"
#include "parameters.h"

/** sequence_length / (this number) is a measure of how hard the
    alignment code will work to find the optimal alignment; in fact
    this gives a worst case bound on the number of loop iterations */
#define GREEDY_MAX_COST_FRACTION 2

/** The largest distance to be examined for an optimal alignment */
#define GREEDY_MAX_COST 2000

/** Bookkeeping structure for greedy alignment. When aligning
    two sequences, the members of this structure store the
    largest offset into the second sequence that leads to a
    high-scoring alignment for a given start point */
struct SGreedyOffset {
    Int4 insert_off;    /**< Best offset for a path ending in an insertion */
    Int4 match_off;     /**< Best offset for a path ending in a match */
    Int4 delete_off;    /**< Best offset for a path ending in a deletion */
};

struct SGreedySeed {
    Int4 start_q;       /**< query offset of start of run of matches */
    Int4 start_s;       /**< subject offset of start of run of matches */
    Int4 match_length;  /**< length of run of matches */
};

/** Space structure for greedy alignment algorithm */
struct SMBSpace {
    SGreedyOffset* space_array; /**< array of bookkeeping structures */
    Int4 space_allocated;       /**< number of structures allocated */
    Int4 space_used;            /**< number of structures actually in use */
    struct SMBSpace *next;      /**< pointer to next structure in list */
};

class MBSpaceMgr
{
public:
  SMBSpace* space;
  
public:
  SMBSpace* MBSPaceNew(int num_space_arrays);
  void RefreshMBSpace();
  void MBSpaceFree();
  SGreedyOffset* GetMBSpace(Int4 num_alloc);
  MBSpaceMgr();
  ~MBSpaceMgr();
};

/** All auxiliary memory needed for the greedy extension algorithm. */
struct SGreedyAlignMem {
   Int4 max_dist;           /* <max distance to search */
   Int4 xdrop;              /* Xdrop value */
   Int4** last_seq2_off;              /**< 2-D array of distances */
   Int4* max_score;                   /**< array of maximum scores */
   SGreedyOffset** last_seq2_off_affine;  /**< Like last_seq2_off but for
                                               affine searches */
   Int4* diag_bounds;                 /**< bounds on ranges of diagonals */
   //SMBSpace* space;                   /**< local memory pool for SGreedyOffset structs */
   MBSpaceMgr* spacemgr;
   
   SGreedyAlignMem() : last_seq2_off(NULL), max_score(NULL), last_seq2_off_affine(NULL), diag_bounds(NULL), spacemgr(NULL) {}
   void GreedyAlignMemFree();
   int MemoryAlloc(const BlastScoringParameters* score_params,
                   const BlastExtensionParameters* ext_params,
                   Int4 arg_max_d, Int4 arg_Xdrop);
   ~SGreedyAlignMem();
};

#endif	/* GREEDY_ALIGN_H */

