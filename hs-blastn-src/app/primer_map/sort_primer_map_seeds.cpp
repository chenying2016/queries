#include "sort_primer_map_seeds.h"

#include <boost/algorithm/sorting/spread_sort.hpp>

using namespace boost;

/////////////////

struct PrimerMapWord_HashLessThan {
    inline bool operator()(const PrimerMapWord& x, const PrimerMapWord& y) const {
        return x.hash < y.hash;
    }
};

struct PrimerMapWord_HashRightShift {
    inline u32 operator()(const PrimerMapWord& x, const unsigned offset) const {
        return x.hash >> offset;
    }
};

extern "C"
void
sort_pm_word_hash_lt(size_t n, PrimerMapWord* a)
{
    integer_sort(a, a + n, PrimerMapWord_HashRightShift(), PrimerMapWord_HashLessThan());
}

///////////////////

struct PrimerMapWord_KeyLessThan {
    inline bool operator()(const PrimerMapWord& x, const PrimerMapWord& y) const {
        return x.key < y.key;
    }
};

struct PrimerMapWord_KeyRightShift {
    inline u64 operator()(const PrimerMapWord& x, const unsigned offset) const {
        return x.key >> offset;
    }
};

extern "C"
void
sort_pm_word_key_lt(size_t n, PrimerMapWord* a)
{
    integer_sort(a, a + n, PrimerMapWord_KeyRightShift(), PrimerMapWord_KeyLessThan());
}

//////////////////

struct PrimerMapSeed_SoffLessThan {
    inline bool operator()(const PrimerMapSeed& x, const PrimerMapSeed& y) const {
        return x.subject_offset < y.subject_offset;
    }
};

struct PrimerMapSeed_SoffRightShift {
    inline int operator()(const PrimerMapSeed& x, const unsigned offset) const {
        return x.subject_offset >> offset;
    }
};

extern "C"
void
sort_pm_seed_soff_lt(size_t n, PrimerMapSeed* a)
{
    integer_sort(a, a + n, PrimerMapSeed_SoffRightShift(), PrimerMapSeed_SoffLessThan());
}

/////////////////

struct PrimerMapSeed_KeyLessThan {
    inline bool operator()(const PrimerMapSeed& x, const PrimerMapSeed& y) const {
        return x.key < y.key;
    }
};

struct PrimerMapSeed_KeyRightShift {
    inline u64 operator()(const PrimerMapSeed& x, const unsigned offset) const {
        return x.key >> offset;
    }
};

extern "C"
void
sort_pm_seed_key_lt(size_t n, PrimerMapSeed* a)
{
    integer_sort(a, a + n, PrimerMapSeed_KeyRightShift(), PrimerMapSeed_KeyLessThan());
}

/////////////////