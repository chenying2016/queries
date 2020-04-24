/* $Id: greedy_align.c 392014 2013-03-13 14:36:38Z maning $
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
 * Functions to perform greedy affine and non-affine gapped alignment
 *
 */

/** @file greedy_align.c
 * Functions to perform greedy affine and non-affine gapped alignment.
 * Reference:
 * Zhang et. al., "A Greedy Algorithm for Aligning DNA Sequences"
 * Journal of Computational Biology vol 7 pp 203-214
 */

#include "greedy_align.h"
#include "math.h"

#include <cstdlib>
#include <cassert>
#include <cstdio>

SMBSpace* MBSpaceMgr::MBSPaceNew(int num_space_arrays)
{
  SMBSpace* new_space;
  const Int4 kMinSpace = 1000000;
  
  num_space_arrays = MAX(kMinSpace, num_space_arrays);
  
  new_space = (SMBSpace*)malloc(sizeof(SMBSpace));
  assert(new_space != NULL);
  
  new_space->space_array = (SGreedyOffset*)malloc(num_space_arrays * sizeof(SGreedyOffset));
  assert(new_space->space_array != NULL);
  
  new_space->space_used = 0;
  new_space->space_allocated = num_space_arrays;
  new_space->next = NULL;
  
  if (0)
  {
      fprintf(stderr, "[%s] allocate memory: %lu\n", __func__,
              sizeof(SGreedyOffset) * num_space_arrays);
  }
  
  return new_space;
}

void MBSpaceMgr::RefreshMBSpace()
{
  SMBSpace* p = space;
  while(p != NULL)
  {
    p->space_used = 0;
    p = p->next;
  } 
}

void MBSpaceMgr::MBSpaceFree()
{
  SMBSpace* p = space;
  SMBSpace* next_space;
  while(p != NULL)
  {
    next_space = p->next;
    free(p->space_array);
    free(p);
    p = next_space;
  }
  space = NULL;
}

MBSpaceMgr::MBSpaceMgr()
{
  space = NULL;
}

MBSpaceMgr::~MBSpaceMgr()
{
  MBSpaceFree();
}

SGreedyOffset* MBSpaceMgr::GetMBSpace(Int4 num_alloc)
{
  SGreedyOffset* out_ptr;
  if (num_alloc < 0)
    return NULL;
  
  SMBSpace* p = space;
  
  while (p->space_used + num_alloc > p->space_allocated) {
    if (p->next == NULL) {
      p->next = MBSPaceNew(num_alloc);
      if (p->next == NULL) return NULL;
    }
    p = p->next;
  }
  
  out_ptr = p->space_array + p->space_used;
  p->space_used += num_alloc;
  
  return out_ptr;  
}

void SGreedyAlignMem::GreedyAlignMemFree()
{
  if (last_seq2_off)
  {
    free(last_seq2_off[0]); last_seq2_off[0] = NULL;
    free(last_seq2_off); last_seq2_off = NULL;
  }
  else
  {
    if (last_seq2_off_affine) {
      free(last_seq2_off_affine[0]); last_seq2_off_affine[0] = NULL;
      free(last_seq2_off_affine); last_seq2_off_affine = NULL;
    }
    free(diag_bounds); diag_bounds = NULL;
  }
  free(max_score); max_score = NULL;
  if (spacemgr) delete spacemgr;
}

int SGreedyAlignMem::MemoryAlloc(const BlastScoringParameters* score_params,
                                 const BlastExtensionParameters* ext_params,
                                 Int4 arg_max_d, Int4 arg_Xdrop)
{
  GreedyAlignMemFree();  
    
  Int4 max_d_1, d_diff, max_cost, gd, i;
  Int4 reward, penalty, gap_open, gap_extend;
  Int4 Mis_cost, GE_cost;

  if (score_params == NULL || (!ext_params && !arg_Xdrop))
    return 1;

    if (score_params->reward % 2 == 1) {
      reward = 2 * score_params->reward;
      penalty = -2 * score_params->penalty;
      if (!arg_Xdrop)
        arg_Xdrop = 2 * MAX(ext_params->gap_x_dropoff, ext_params->gap_x_dropoff_final);
      gap_open = 2 * score_params->gap_open;
      gap_extend = 2 * score_params->gap_extend;
    } else {
      reward = score_params->reward;
      penalty = score_params->penalty;
      if (!arg_Xdrop) {
        arg_Xdrop = MAX(ext_params->gap_x_dropoff, ext_params->gap_x_dropoff_final);
      }
      gap_open = score_params->gap_open;
      gap_extend = score_params->gap_extend;
    }

    if (gap_open == 0 && gap_extend == 0)
    {
      gap_extend = reward / 2 + penalty;
    }

    max_dist = arg_max_d;
    xdrop = arg_Xdrop;

    if (score_params->gap_open == 0 && score_params->gap_extend == 0) {
      d_diff = (xdrop + reward / 2) / (penalty + reward) + 1;

      last_seq2_off = (Int4**)malloc((arg_max_d + 2) * sizeof(Int4*));
      if (last_seq2_off == NULL)
        return 1;
      last_seq2_off[0] = (Int4*)malloc((arg_max_d + arg_max_d + 6) * sizeof(Int4) * 2);
      if (last_seq2_off[0] == NULL)
        return 1;

      last_seq2_off[1] = last_seq2_off[0] + max_dist + max_dist + 6;
      last_seq2_off_affine = NULL;
      diag_bounds = NULL;
    } else
    {
      last_seq2_off = NULL;
      Mis_cost = reward + penalty;
      GE_cost = gap_extend + reward / 2;
      max_d_1 = arg_max_d;
      arg_max_d *= GE_cost;
      max_cost = MAX(Mis_cost, gap_open + GE_cost);
      gd = Gdb3(&Mis_cost, &gap_open, &GE_cost);
      d_diff = (xdrop + reward / 2)/gd + 1;
      diag_bounds = (Int4*)calloc(2 * (arg_max_d + 1 + max_cost), sizeof(Int4));
      last_seq2_off_affine = (SGreedyOffset**)malloc((MAX(arg_max_d, max_cost) + 2) * sizeof(SGreedyOffset*));
      if (!diag_bounds || !last_seq2_off_affine) return 1;

      last_seq2_off_affine[0] = (SGreedyOffset*)calloc((2 * max_d_1 + 6), sizeof(SGreedyOffset) * (max_cost + 1));
      for (i = 1; i <= max_cost; ++i)
      {
        last_seq2_off_affine[i] = last_seq2_off_affine[i-1] + 2 * max_d_1 + 6;
      }

      if (!last_seq2_off_affine || !last_seq2_off_affine[0]) return 1;
    }

    max_score = (Int4*)malloc(sizeof(Int4) * (arg_max_d + 1 + d_diff));

    spacemgr = new MBSpaceMgr();
    spacemgr->space = spacemgr->MBSPaceNew(0);   

    return 0;
}

SGreedyAlignMem::~SGreedyAlignMem()
{
  GreedyAlignMemFree();
}
