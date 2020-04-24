#include "subseq_hit.h"

#include <algorithm>

using namespace std;

extern "C"
void sort_subseq_hit_score_gt(size_t n, HbnSubseqHit* a)
{
    sort(a, 
        a + n, 
        [](const HbnSubseqHit& lhs, const HbnSubseqHit& rhs){ return lhs.score > rhs.score; } );
}

extern "C"
void sort_subseq_hit_sid_lt(size_t n, HbnSubseqHit* a)
{
    sort(a,
        a + n,
        [](const HbnSubseqHit& lhs, const HbnSubseqHit& rhs){ return lhs.sid < rhs.sid; });
}

extern "C"
void sort_subseq_hit_sfrom_lt(size_t n, HbnSubseqHit* a)
{
    sort(a,
        a + n,
        [](const HbnSubseqHit& lhs, const HbnSubseqHit& rhs){ return (lhs.sfrom < rhs.sfrom) || (lhs.sfrom == rhs.sfrom && lhs.sto > rhs.sto); });
}