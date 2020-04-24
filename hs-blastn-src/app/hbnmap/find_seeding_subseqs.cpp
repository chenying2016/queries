#include "find_seeding_subseqs.h"

#include "symdust.hpp"

extern "C"
int find_seeding_subseqs(const unsigned char* seq, 
        const unsigned int seq_size,
        const int min_size,
        vec_int_pair* good_regions)
{
    CSymDustMasker masker;
    std::vector< std::pair<CSymDustMasker::size_type, CSymDustMasker::size_type> > masked_locs;
    masker.GetMaskedLocs(seq, seq_size, masked_locs);
    IntPair ip;
    kv_clear(*good_regions);
    if (masked_locs.empty()) {
        ip.first = 0;
        ip.second = seq_size;
        kv_push(IntPair, *good_regions, ip);
        return 1;
    }

    CSymDustMasker::size_type from = 0, to;
    for (size_t i = 0; i < masked_locs.size(); ++i) {
        to = masked_locs[i].first;
        if (to - from >= min_size) {
            ip.first = from;
            ip.second = to;
            kv_push(IntPair, *good_regions, ip);
        }
        from = masked_locs[i].second + 1;
    }

    if (seq_size - from >= min_size) {
        ip.first = from;
        ip.second = seq_size;
        kv_push(IntPair, *good_regions, ip);
    }

#if 0
    printf("Bad Regions:\n");
    for (size_t i = 0; i < masked_locs.size(); ++i) {
        printf("[%d, %d]\n", masked_locs[i].first, masked_locs[i].second);
    }
    printf("Good Regions:\n");
    for (size_t i = 0; i < kv_size(*good_regions); ++i) {
        printf("[%d, %d]\n", kv_A(*good_regions, i).first, kv_A(*good_regions, i).second);
    }
#endif

    return kv_size(*good_regions);
}