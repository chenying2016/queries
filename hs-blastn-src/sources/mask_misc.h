#ifndef MASK_MISC_H
#define MASK_MISC_H

#include <seq_masker.hpp>
#include <algorithm>

void SeqLocReverse(CSeqMasker::TMaskList& masks, Int4 query_length);

void SeqLocCombine(CSeqMasker::TMaskList& mask_loc, Int4 link_value);

void ComplementMaskLocations(Int4 query_length,
							 Int4 length_threshold,
                             CSeqMasker::TMaskList& masked_locs,
                             CSeqMasker::TMaskList& unmasked_locs);

#endif // MASK_MISC_H
