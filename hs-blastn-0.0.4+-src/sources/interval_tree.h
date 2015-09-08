/* $Id: blast_itree.h 240628 2011-02-09 14:37:10Z coulouri $
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
 * Author: Jason Papadopoulos
 *
 */

/** @file blast_itree.h
 * Interface for an interval tree, used for fast HSP containment tests
 */

/**
 * This file a modified version of blast_itree.h
*/


#ifndef INTERVAL_TREE_H
#define	INTERVAL_TREE_H

#include "gapalign.h"

typedef Int8 ITree_idxt;
typedef HSP* ITree_keyt;
#define ITreeInvalidKey NULL

enum EIntervalDirection
{
    eIntervalTreeLeft,
    eIntervalTreeRight,
    eIntervalTreeNeither
};

struct SIntervalNode
{
    ITree_idxt leftend;
    ITree_idxt rightend;
    Int4 leftptr;
    Int4 rightptr;
    Int4 midptr;
    ITree_keyt key;
};

class IntervalTree
{
private:
    SIntervalNode* nodes;
    Int4 num_alloc;
    Int4 num_used;
    ITree_idxt s_min;
    ITree_idxt s_max;

private:
    Int4 IntervalNodeInit(Int4 parent_index, enum EIntervalDirection dir, Int2& retval);
    Int4 IntervalRootNodeInit(ITree_idxt region_start, ITree_idxt region_end, Int2& retval);

private:
    const HSP* HSPsHaveCommonEndpoint(const HSP* in_hsp, const HSP* tree_hsp, enum EIntervalDirection which_end);
    Boolean MidpointTreeHasHSPEndpoint(Int4 root_index, const HSP* in_hsp, enum EIntervalDirection which_end);
    Boolean IntervalTreeHasHSPEndpoint(const HSP* in_hsp, enum EIntervalDirection which_end);
    Boolean HSPIsContained(const HSP* in_hsp, const HSP* tree_hsp, Int4 min_diag_seperation);
    Boolean MidpointTreeContainsHSP(Int4 root_index, const HSP* in_hsp, Int4 min_diag_seperation);

public:
    IntervalTree(ITree_idxt qs = 0, ITree_idxt qe = 0, ITree_idxt ss = 0, ITree_idxt se = 0);
    void Destroy();
    ~IntervalTree();

public:
    void Reset(ITree_idxt qs, ITree_idxt qe, ITree_idxt ss, ITree_idxt se);
    Boolean IntervalTreeContainsHSP(const HSP* in_hsp, Int4 min_diag_seperation);
    Int2 IntervalTreeAddHSP(HSP* hsp);
};

#endif	/* INTERVAL_TREE_H */

