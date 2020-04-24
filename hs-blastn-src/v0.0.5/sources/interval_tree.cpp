/* $Id: blast_itree.c 358499 2012-04-03 14:48:04Z coulouri $
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

/** @file blast_itree.c
 * Functions that implement an interval tree for fast HSP containment tests
 */

/**
 * This file is a modified version of blast_itree.c.
*/

#include "interval_tree.h"

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>

Int4 IntervalTree::IntervalNodeInit(Int4 parent_index, enum EIntervalDirection dir, Int2& ret_status)
{
    Int4 new_index;
    ITree_idxt midpt;
    SIntervalNode *new_node;
    SIntervalNode *parent_node;

    ret_status = 0;

    if (num_used == num_alloc) {
        num_alloc = 2 * num_alloc;
        nodes = (SIntervalNode *)realloc(nodes, num_alloc * sizeof(SIntervalNode));
    }

    if(nodes == NULL)
    {
         ret_status = 1;
         return 0;
    }

    new_index = num_used++;
    new_node = nodes + new_index;
    new_node->key = ITreeInvalidKey;

    if (dir == eIntervalTreeNeither)
        return new_index;

    /* fields in the node are only filled in if a parent
       node is specified */

    parent_node = nodes + parent_index;
    new_node->leftptr = 0;
    new_node->midptr = 0;
    new_node->rightptr = 0;
    midpt = (parent_node->leftend + parent_node->rightend) / 2;

    /* the endpoints of the new node depend on whether
       it is for the left or right subtree of the parent.
       The two subregions do not overlap, may be of length
       one, and must completely cover the parent region */

    if (dir == eIntervalTreeLeft) {
        new_node->leftend = parent_node->leftend;
        new_node->rightend = midpt;
    }
    else {
        new_node->leftend = midpt + 1;
        new_node->rightend = parent_node->rightend;
    }

    return new_index;
}

Int4 IntervalTree::IntervalRootNodeInit(ITree_idxt region_start, ITree_idxt region_end, Int2& retval)
{
    Int4 new_index;
    SIntervalNode *new_node;

    new_index = IntervalNodeInit(0, eIntervalTreeNeither, retval);
    if(retval != 0)
       return 0;

    new_node = nodes + new_index;
    new_node->leftptr = 0;
    new_node->midptr = 0;
    new_node->rightptr = 0;
    new_node->key = ITreeInvalidKey;
    new_node->leftend = region_start;
    new_node->rightend = region_end;
    return new_index;
}

void IntervalTree::Destroy()
{
    if (nodes)
    {
        free(nodes);
        nodes = NULL;
        num_alloc = num_used = 0;
    }
}

IntervalTree::~IntervalTree()
{
    Destroy();
}

void IntervalTree::Reset(ITree_idxt qs, ITree_idxt qe, ITree_idxt ss, ITree_idxt se)
{
    num_used = 0;
    s_min = ss;
    s_max = se;
    Int2 retval;
    IntervalRootNodeInit(qs, qe, retval);
    assert(retval == 0);
}

IntervalTree::IntervalTree(ITree_idxt qs, ITree_idxt qe, ITree_idxt ss, ITree_idxt se)
{
    const Int4 kSize = 256;
    Int2 retval = 0;

    nodes = (SIntervalNode*)calloc(kSize, sizeof(SIntervalNode));
    assert(nodes != NULL);

    num_alloc = kSize;
    num_used = 0;
    s_min = ss;
    s_max = se;

    IntervalRootNodeInit(qs, qe, retval);
    assert(retval == 0);
}

/** TRUE if c is between a and b; f between d and e.  Determines if the
 * coordinates are already in an HSP that has been evaluated.
*/
#define CONTAINED_IN_HSP(a,b,c,d,e,f) \
    (((a <= c && b >= c) && (d <= f && e >= f)) ? TRUE : FALSE)

/** Are the two HSPs within a given number of diagonals from each other? */
#define MB_HSP_CLOSE(q1, s1, q2, s2, c) \
(ABS(((q1)-(s1)) - ((q2)-(s2))) < c)

Boolean IntervalTree::HSPIsContained(const HSP* in_hsp, const HSP* tree_hsp, Int4 min_diag_seperation)
{
    if (in_hsp->subject_id != tree_hsp->subject_id) return FALSE;

    if (in_hsp->score <= tree_hsp->score &&
        CONTAINED_IN_HSP(tree_hsp->q_off, tree_hsp->q_end, in_hsp->q_off,
                         tree_hsp->s_off, tree_hsp->s_end, in_hsp->s_off) &&
        CONTAINED_IN_HSP(tree_hsp->q_off, tree_hsp->q_end, in_hsp->q_end,
                         tree_hsp->s_off, tree_hsp->s_end, in_hsp->s_end))
    {
        if (min_diag_seperation == 0)
        return TRUE;        

        if (MB_HSP_CLOSE(tree_hsp->q_off, tree_hsp->s_off,
                         in_hsp->q_off, in_hsp->s_off,
                         min_diag_seperation) ||
            MB_HSP_CLOSE(tree_hsp->q_end, tree_hsp->s_end,
                         in_hsp->q_end, in_hsp->s_end,
                         min_diag_seperation))
        {
            return TRUE;
        }
    }
    return FALSE;
}

const HSP* IntervalTree::HSPsHaveCommonEndpoint(const HSP* in_hsp, const HSP* tree_hsp,
                                                enum EIntervalDirection which_end)
{
    Boolean match;

    if (which_end == eIntervalTreeLeft) {
        match = in_hsp->q_off== tree_hsp->q_off &&
                in_hsp->s_off == tree_hsp->s_off;
    }
    else {
        match = in_hsp->q_end == tree_hsp->q_end &&
                in_hsp->s_end == tree_hsp->s_end;
    }

    if (match) {
        ITree_idxt in_q_length, tree_q_length, in_s_length, tree_s_length;

        /* keep the higher scoring HSP */

        if (in_hsp->score > tree_hsp->score)
            return in_hsp;
        if (in_hsp->score < tree_hsp->score)
            return tree_hsp;

        /* for equal scores, pick the shorter HSP */
        in_q_length = in_hsp->q_end - in_hsp->q_off;
        tree_q_length = tree_hsp->q_end - tree_hsp->q_off;
        if (in_q_length > tree_q_length)
            return tree_hsp;
        if (in_q_length < tree_q_length)
            return in_hsp;

        in_s_length = in_hsp->s_end - in_hsp->s_off;
        tree_s_length = tree_hsp->s_end - tree_hsp->s_off;
        if (in_s_length > tree_s_length)
            return tree_hsp;
        if (in_s_length < tree_s_length)
            return in_hsp;

        /* HSPs are identical; favor the one already in the tree */

        return tree_hsp;
    }

    return NULL;
}


Boolean
IntervalTree::MidpointTreeHasHSPEndpoint(Int4 root_index, const HSP* in_hsp,
                                         enum EIntervalDirection which_end)
{
    SIntervalNode *root_node = nodes + root_index;
    SIntervalNode *list_node, *next_node;
    Int4 tmp_index;
    ITree_idxt target_offset;
    ITree_idxt midpt;

    if (which_end == eIntervalTreeLeft)
        target_offset = in_hsp->s_off;
    else
        target_offset = in_hsp->s_end;

    /* Descend the tree */

    while (1) {

        ASSERT(target_offset >= root_node->leftend);
        ASSERT(target_offset <= root_node->rightend);

        /* First perform matching endpoint tests on all of the HSPs
           in the midpoint list for the current node. If the input
           shares an endpoint with an HSP already in the list, and the
           HSP in the list is 'better', signal that in_hsp should not
           be added to the tree later. Otherwise remove matching HSPs
           from the list. */

        tmp_index = root_node->midptr;
        list_node = root_node;
        next_node = nodes + tmp_index;
        while (tmp_index != 0) {
        	HSP* tree_hsp = next_node->key;
            const HSP *best_hsp = HSPsHaveCommonEndpoint(in_hsp,
            						tree_hsp,
                                                        which_end);

            tmp_index = next_node->midptr;
            if (best_hsp == tree_hsp)
                return TRUE;
            else if (best_hsp == in_hsp)
                list_node->midptr = tmp_index;

            list_node = next_node;
            next_node = nodes + tmp_index;
        }

        /* Descend to the left or right subtree, whichever one
           contains the endpoint from in_hsp */

        tmp_index = 0;
        midpt = (root_node->leftend + root_node->rightend) / 2;
        if (target_offset < midpt)
            tmp_index = root_node->leftptr;
        else if (target_offset > midpt)
            tmp_index = root_node->rightptr;

        /* If there is no such subtree, then all of the HSPs that
           could possibly have a common endpoint with it have already
           been examined */

        if (tmp_index == 0)
            return FALSE;

        next_node = nodes + tmp_index;
        
        if (next_node->key != ITreeInvalidKey) {

            /* reached a leaf; compare in_hsp with the alignment
               in the leaf. Whether or not there's a match, traversal
               is finished */

        	HSP* tree_hsp = next_node->key;
            const HSP *best_hsp = HSPsHaveCommonEndpoint(in_hsp,
                                                 tree_hsp,
                                                 which_end);
            if (best_hsp == tree_hsp) {
                return TRUE;
            }
            else if (best_hsp == in_hsp) {
                /* leaf gets removed */
                if (target_offset < midpt)
                    root_node->leftptr = 0;
                else if (target_offset > midpt)
                    root_node->rightptr = 0;
                return FALSE;
            }
            break;
        }
        root_node = next_node;          /* descend to next node */
    }
    return FALSE;
}

Boolean IntervalTree::IntervalTreeHasHSPEndpoint(const HSP* in_hsp, enum EIntervalDirection which_end)
{
    SIntervalNode *root_node = nodes;
    SIntervalNode *next_node;
    Int4 tmp_index;
    ITree_idxt target_offset;
    ITree_idxt midpt;

    if (which_end == eIntervalTreeLeft)
        target_offset = in_hsp->q_off;
    else
        target_offset = in_hsp->q_end;

    /* Descend the tree */

    while (1) {

        ASSERT(target_offset >= root_node->leftend);
        ASSERT(target_offset <= root_node->rightend);

        /* First perform matching endpoint tests on all of the HSPs
           in the midpoint tree for the current node */

        tmp_index = root_node->midptr;
        if (tmp_index != 0) {
            if (MidpointTreeHasHSPEndpoint(tmp_index, in_hsp,
                                            which_end)) {
                return TRUE;
            }
        }

        /* Descend to the left or right subtree, whichever one
           contains the endpoint from in_hsp */

        tmp_index = 0;
        midpt = (root_node->leftend + root_node->rightend) / 2;
        if (target_offset < midpt)
            tmp_index = root_node->leftptr;
        else if (target_offset > midpt)
            tmp_index = root_node->rightptr;

        /* If there is no such subtree, or the HSP straddles the center
           of the current node, then all of the HSPs that could possibly
           have a common endpoint with it have already been examined */

        if (tmp_index == 0)
            return FALSE;

        next_node = nodes + tmp_index;
        if (next_node->key != ITreeInvalidKey) {

            /* reached a leaf; compare in_hsp with the alignment
               in the leaf. Whether or not there's a match, traversal
               is finished */

        	HSP* tree_hsp = next_node->key;
            const HSP* best_hsp = HSPsHaveCommonEndpoint(in_hsp,
                                              tree_hsp,
                                              which_end);
            if (best_hsp == tree_hsp) {
                return TRUE;
            }
            else if (best_hsp == in_hsp) {
                /* leaf gets removed */
                if (target_offset < midpt)
                    root_node->leftptr = 0;
                else if (target_offset > midpt)
                    root_node->rightptr = 0;
                return FALSE;
            }
            break;
        }
        root_node = next_node;          /* descend to next node */
    }
    return FALSE;
}

Boolean IntervalTree::MidpointTreeContainsHSP(Int4 root_index, const HSP* in_hsp, Int4 min_diag_separation)
{
    SIntervalNode *node = nodes + root_index;
    ITree_idxt region_start = in_hsp->s_off;
    ITree_idxt region_end = in_hsp->s_end;
    ITree_idxt middle;
    Int4 tmp_index = 0;

    /* Descend the tree */

    while (node->key == ITreeInvalidKey) {

        ASSERT(region_start >= node->leftend);
        ASSERT(region_end <= node->rightend);

        /* First perform containment tests on all of the HSPs
           in the midpoint list for the current node. These
           HSPs are not indexed in a tree format, so all HSPs
           in the list must be examined */

        tmp_index = node->midptr;
        while (tmp_index != 0) {
            SIntervalNode *tmp_node = nodes + tmp_index;

            HSP* tree_hsp = tmp_node->key;
            if (HSPIsContained(in_hsp,
                                 tree_hsp,
                                 min_diag_separation)) {
                return TRUE;
            }
            tmp_index = tmp_node->midptr;
        }

        /* Descend to the left subtree if the input HSP lies completely
           to the left of this node's center, or to the right subtree if
           it lies completely to the right */

        tmp_index = 0;
        middle = (node->leftend + node->rightend) / 2;
        if (region_end < middle)
            tmp_index = node->leftptr;
        else if (region_start > middle)
            tmp_index = node->rightptr;

        /* If there is no such subtree, or the HSP straddles the center
           of the current node, then all of the HSPs that could possibly
           contain it have already been examined */

        if (tmp_index == 0)
            return FALSE;

        node = nodes + tmp_index;
    }

    /* Reached a leaf of the tree */

    return HSPIsContained(in_hsp,
                          node->key,
                          min_diag_separation);
}

Boolean IntervalTree::IntervalTreeContainsHSP(const HSP* hsp, Int4 min_diag_separation)
{
    SIntervalNode *node = nodes;
    ITree_idxt query_start = 0; //s_GetQueryStrandOffset(query_info, hsp->context);
    ITree_idxt region_start = query_start + hsp->q_off;
    ITree_idxt region_end = query_start + hsp->q_end;
    ITree_idxt middle;
    Int4 tmp_index = 0;

    ASSERT(region_start >= node->leftend);
    ASSERT(region_end <= node->rightend);
    ASSERT(hsp->s_off >= s_min);
    ASSERT(hsp->s_end <= s_max);
    ASSERT(hsp->q_off <= hsp->q_end);
    ASSERT(hsp->s_off <= hsp->s_end);

    /* Descend the tree */

    while (node->key == ITreeInvalidKey) {

        ASSERT(region_start >= node->leftend);
        ASSERT(region_end <= node->rightend);

        /* First perform containment tests on all of the HSPs
           in the midpoint tree for the current node */

        tmp_index = node->midptr;
        if (tmp_index > 0) {
            if (MidpointTreeContainsHSP(tmp_index,
                                        hsp,
                                        min_diag_separation)) {
                return TRUE;
            }
        }

        /* Descend to the left subtree if the input HSP lies completely
           to the left of this node's center, or to the right subtree if
           it lies completely to the right */

        tmp_index = 0;
        middle = (node->leftend + node->rightend) / 2;
        if (region_end < middle)
            tmp_index = node->leftptr;
        else if (region_start > middle)
            tmp_index = node->rightptr;

        /* If there is no such subtree, or the HSP straddles the center
           of the current node, then all of the HSPs that could possibly
           contain it have already been examined */

        if (tmp_index == 0)
            return FALSE;

        node = nodes + tmp_index;
    }

    /* Reached a leaf of the tree */

    return HSPIsContained(hsp,
                          node->key,
                          min_diag_separation);
}

Int2 IntervalTree::IntervalTreeAddHSP(HSP* hsp)
{
    ITree_idxt query_start;
    ITree_idxt old_region_start;
    ITree_idxt old_region_end;
    ITree_idxt region_start;
    ITree_idxt region_end;
    HSP *old_hsp;
    Int4 root_index;
    Int4 new_index;
    Int4 mid_index;
    Int4 old_index;
    ITree_idxt middle;
    enum EIntervalDirection which_half;
    Boolean index_subject_range = FALSE;
    Int2 retval = 0;
    Int4 q_start;
    Int4 mid_index2;

    /* Determine the query strand containing the input HSP.
       Only the strand matters for containment purposes,
       not the precise value of the query frame */

    query_start = 0; //s_GetQueryStrandOffset(query_info, hsp->context);

	region_start = query_start + hsp->q_off;
	region_end = query_start + hsp->q_end;

    ASSERT(region_start >= nodes->leftend);
    ASSERT(region_end <= nodes->rightend);
    ASSERT(hsp->q_off <= hsp->q_end);
    ASSERT(hsp->s_off <= hsp->s_end);

	ASSERT(hsp->s_off >= s_min);
	ASSERT(hsp->s_end <= s_max);

	/* Before adding the HSP, determine whether one or more
	 HSPs already in the tree share a common endpoint with
	 in_hsp. Remove from the tree any leaves containing
	 such an HSP whose score is lower than in_hsp.

	 Note that in_hsp might share an endpoint with a
	 higher-scoring HSP already in the tree, in which case
	 in_hsp should not be added. There is thus a possibility
	 that in_hsp will remove an alignment from the tree and
	 then another alignment will remove in_hsp. This is arguably
	 not the right behavior, but since the tree is only for
	 containment tests the worst that can happen is that
	 a rare extra gapped alignment will be computed */

	if (IntervalTreeHasHSPEndpoint(hsp, eIntervalTreeLeft))
	{
		return retval;
	}
	if (IntervalTreeHasHSPEndpoint(hsp, eIntervalTreeRight))
	{
		return retval;
	}


    /* begin by indexing the HSP query offsets */

    index_subject_range = FALSE;

    /* encapsulate the input HSP in an SIntervalNode */
    root_index = 0;
    new_index = IntervalNodeInit(0, eIntervalTreeNeither, retval);
    if (retval)
         return retval;
    nodes[new_index].leftptr = query_start;
    nodes[new_index].midptr = 0;
    nodes[new_index].key = hsp;

    /* Descend the tree to reach the correct subtree for the new node */

    while (1) {

        ASSERT(region_start >= nodes[root_index].leftend);
        ASSERT(region_end <= nodes[root_index].rightend);

        middle = (nodes[root_index].leftend +
                  nodes[root_index].rightend) / 2;

        if (region_end < middle) {

            /* new interval belongs in left subtree. If there
               are no leaves in that subtree, finish up */

            if (nodes[root_index].leftptr == 0) {
                nodes[root_index].leftptr = new_index;
                return retval;
            }

            /* A node is already in this subtree. If it is not a
               leaf node, descend to it and analyze in the next
               loop iteration. Otherwise, schedule the subtree to
               be split */

            old_index = nodes[root_index].leftptr;
            if (nodes[old_index].key == ITreeInvalidKey) {
                root_index = old_index;
                continue;
            }
            else {
                which_half = eIntervalTreeLeft;
            }
        }
        else if (region_start > middle) {

            /* new interval belongs in right subtree. If there
               are no leaves in that subtree, finish up */

            if (nodes[root_index].rightptr == 0) {
                nodes[root_index].rightptr = new_index;
                return retval;
            }

            /* A node is already in this subtree. If it is not a
               leaf node, descend to it and analyze in the next
               loop iteration. Otherwise, schedule the subtree to
               be split */

            old_index = nodes[root_index].rightptr;
            if (nodes[old_index].key == ITreeInvalidKey) {
                root_index = old_index;
                continue;
            }
            else {
                which_half = eIntervalTreeRight;
            }
        }
        else {

            /* the new interval crosses the center of the node, and
               so has a "shadow" in both subtrees */

            // Added support for eQueryOnlyStrandIndifferent -RMH-
            if (index_subject_range) {

                /* If indexing subject offsets already, prepend the
                   new node to the list of "midpoint" nodes and return.
                   midptr is always a linked list if only the query
                   offsets are indexed */

                nodes[new_index].midptr = nodes[root_index].midptr;
                nodes[root_index].midptr = new_index;
                return retval;
            }
            else {

                /* Begin another tree at root_index, that indexes
                   the subject range of the input HSP */

                index_subject_range = TRUE;

                if (nodes[root_index].midptr == 0) {
                    mid_index = IntervalRootNodeInit(s_min,
                                                     s_max, retval);
                    if (retval)
                      return retval;
                    nodes[root_index].midptr = mid_index;
                }
                root_index = nodes[root_index].midptr;

                /* switch from the query range of the input HSP
                   to the subject range */

                region_start = hsp->s_off;
                region_end = hsp->s_end;
                continue;
            }
        }

        /* There are two leaves in the same subtree. Add another
           internal node, reattach the old leaf, and loop

           First allocate the new node. Update the pointer to
           the pool of nodes, since it may change */

        mid_index = IntervalNodeInit(root_index, which_half, retval);
        if (retval)
          return retval;
        old_hsp = nodes[old_index].key;

        /* attach the new internal node */

        if (which_half == eIntervalTreeLeft)
                nodes[root_index].leftptr = mid_index;
        else
                nodes[root_index].rightptr = mid_index;

        /* descend to the new internal node, and attach the old
           leaf to it. The next loop iteration will have to deal
           with attaching the *new* leaf */

        if (index_subject_range) {
            old_region_start = old_hsp->s_off;
            old_region_end = old_hsp->s_end;
        }
        else {
			old_region_start = nodes[old_index].leftptr + old_hsp->q_off;
			old_region_end = nodes[old_index].leftptr + old_hsp->q_end;
        }

        root_index = mid_index;
        middle = (nodes[root_index].leftend +
                  nodes[root_index].rightend) / 2;
        if (old_region_end < middle) {

            /* old leaf belongs in left subtree of new node */
            nodes[mid_index].leftptr = old_index;
        }
        else if (old_region_start > middle) {

            /* old leaf belongs in right subtree of new node */
            nodes[mid_index].rightptr = old_index;
        }
        else {

            /* the old leaf straddles both subtrees. If indexing is
               by subject offset, attach the old leaf to the (empty)
               midpoint list of the new node. If still indexing query
               offsets, then a new tree that indexes subject offsets
               must be allocated from scratch, just to accomodate the
               old leaf */

            // Added support for eQueryOnlyStrandIndifferent -RMH-
            if (index_subject_range) {
                nodes[mid_index].midptr = old_index;
            }
            else {
                mid_index2 = IntervalRootNodeInit(s_min,
                                                  s_max, retval);
                if (retval)
                      return retval;
                old_region_start = old_hsp->s_off;
                old_region_end =  old_hsp->s_end;
                nodes[mid_index].midptr = mid_index2;
                middle = (nodes[mid_index2].leftend +
                          nodes[mid_index2].rightend) / 2;

                if (old_region_end < middle)
                    nodes[mid_index2].leftptr = old_index;
                else if (old_region_start > middle)
                    nodes[mid_index2].rightptr = old_index;
                else
                    nodes[mid_index2].midptr = old_index;
            }
        }
    }
    return retval;
}
