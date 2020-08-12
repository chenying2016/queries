#include "sort_sr_hit_seeds.h"

#include <boost/algorithm/sorting/spread_sort.hpp>

using namespace boost;


/////////////////////

struct ChainSeed_SoffLessThan {
    inline bool operator()(const ChainSeed& x, const ChainSeed& y) const { return x.soff < y.soff; }
};

struct ChainSeed_SoffRightShift {
    inline idx operator()(const ChainSeed &x, const unsigned offset) { return x.soff >> offset; }
};

extern "C"
void
sort_chain_seed_soff_lt(size_t n, ChainSeed* a)
{
    integer_sort(a, a + n, ChainSeed_SoffRightShift(), ChainSeed_SoffLessThan());
}