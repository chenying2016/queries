#include <vector>

#include "mask_misc.h"

void SeqLocReverse(CSeqMasker::TMaskList& masks, Int4 query_length)
{
    int n = masks.size();
    int i;
    for (i = 0; i < n; ++i)
    {
        CSeqMasker::TMaskedInterval& interval = masks[i];
        CSeqMasker::size_type left = query_length - 1 - interval.second;
        CSeqMasker::size_type right = query_length - 1 - interval.first;
        interval.first = left;
        interval.second = right;
    }
}

struct SortMaskedInterval
{
    bool operator()(const CSeqMasker::TMaskedInterval& a, const CSeqMasker::TMaskedInterval& b)
    {
        if (a.first <= b.first) return true;
        return false;
    }
};

void SeqLocCombine(CSeqMasker::TMaskList& mask_loc, Int4 link_value)
{
    std::sort(mask_loc.begin(), mask_loc.end(), SortMaskedInterval());     
    
    /* Merge the overlapping elements */
    int i, curr = 0;
    int n = mask_loc.size();
    for (i = 0; i < n - 1; ++i)
    {
        CSeqMasker::TMaskedInterval& next_loc = mask_loc[i+1];
        const Int4 stop = mask_loc[curr].second;

        if ((stop + link_value) > static_cast<int>(next_loc.first))
        {
            mask_loc[curr].second = MAX(stop, static_cast<int>(next_loc.second));
            next_loc.second = 0;
        }
        else
        {
            curr = i + 1;
        }
    }

    /* Rebuild the linked list */
    curr = 1;
    for (i = 1; i < n; ++i)
    {
        if (mask_loc[i].second != 0)
        {
            mask_loc[curr].first = mask_loc[i].first;
            mask_loc[curr].second = mask_loc[i].second;
            ++curr;
        }
    }

    mask_loc.resize(curr);
}

void ComplementMaskLocations(Int4 query_length,
                             Int4 length_threshold,
                             CSeqMasker::TMaskList& masked_locs,
                             CSeqMasker::TMaskList& unmasked_locs)
{
    unmasked_locs.clear();

    bool first = true;
    bool last_interval_open = true;
    Int4 start_offset, end_offset, filter_start, filter_end;
    Int4 left = 0, right;

    start_offset = 0;
    end_offset = query_length - 1;
    ASSERT(start_offset <= end_offset);

    CSeqMasker::TMaskList::const_iterator cmiter;
    for (cmiter = masked_locs.begin(); cmiter != masked_locs.end(); ++cmiter)
    {
        const CSeqMasker::TMaskedInterval& seq_range = *cmiter;
        filter_start = start_offset + seq_range.first;
        filter_end = start_offset + seq_range.second;

        if (first)
        {
            last_interval_open = true;
            first = false;

            if (filter_start > start_offset)
                left = start_offset;
            else
            {
                left = filter_end + 1;
                continue;
            }
        }

        right = filter_start - 1;
        if (right - left + 1 >= length_threshold)
        {
            unmasked_locs.push_back(std::make_pair(left, right));
        }

        if (filter_end >= end_offset)
        {
            last_interval_open = false;
            break;
        }
        else
        {
            left = filter_end + 1;
        }
    }

    if (last_interval_open)
    {
        right = end_offset;
        if (right - left + 1 >= length_threshold)
        {
            unmasked_locs.push_back(std::make_pair(left, right));
        }
    }
}

