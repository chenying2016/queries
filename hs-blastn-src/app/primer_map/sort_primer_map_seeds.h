#ifndef __SORT_PRIMER_MAP_SEEDS_H
#define __SORT_PRIMER_MAP_SEEDS_H

#include "primer_map_hit_finder.h"

#ifdef __cplusplus
extern "C" {
#endif

void
sort_pm_word_hash_lt(size_t n, PrimerMapWord* a);

void
sort_pm_word_key_lt(size_t n, PrimerMapWord* a);

void
sort_pm_seed_soff_lt(size_t n, PrimerMapSeed* a);

void
sort_pm_seed_key_lt(size_t n, PrimerMapSeed* a);

#ifdef __cplusplus
}
#endif

#endif // __SORT_PRIMER_MAP_SEEDS_H