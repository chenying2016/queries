//Details for Templated Spreadsort-based implementations and string_sort

//          Copyright Steven J. Ross 2001 - 2009.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org/ for updates, documentation, and revision history
// See http://www.boost.org/libs/algorithm/sorting for library home page.
      
/*
Some improvements suggested by:
Phil Endecott and Frank Gennari
float_mem_cast fix provided by:
Scott McMurray
*/

#ifndef BOOST_SPREAD_SORT_DETAIL_H
#define BOOST_SPREAD_SORT_DETAIL_H
#include <algorithm>
#include <vector>
#include <cstring>
#include <limits>
#include <functional>
#include <boost/static_assert.hpp>
#include <boost/static_warning.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/algorithm/sorting/constants.hpp>
#include <boost/cstdint.hpp>

namespace boost {
  namespace detail {
    //------------------------------------------------- general details

    //This only works on unsigned data types
    template <typename T>
    inline unsigned 
    rough_log_2_size(const T& input) 
    {
      unsigned result = 0;
      //The && is necessary on some compilers to avoid infinite loops
    //it doesn't significantly impair performance
      while((input >> result) && (result < (8*sizeof(T)))) ++result;
      return result;
    }

    //Gets the maximum size which we'll call spread_sort on to control 
    //worst-case performance
    //This is called for a set of bins, instead of bin-by-bin, 
    //to avoid performance overhead
    //This could be replaced by a lookup table of sizeof(Div_type)*8
    //but this isn't done to keep the code general
    template<unsigned log_mean_bin_size, 
      unsigned log_min_split_count, unsigned log_finishing_count>
    inline size_t
    get_max_count(unsigned log_range)
    {
      const size_t typed_one = 1;
      const unsigned min_size = log_mean_bin_size + log_min_split_count;
      //Assuring that constants have valid settings
      BOOST_STATIC_ASSERT(log_min_split_count <= MAX_SPLITS && log_min_split_count > 0);
      BOOST_STATIC_ASSERT(MAX_SPLITS > 1 && MAX_SPLITS < (8 * sizeof(unsigned)));
      BOOST_STATIC_ASSERT(MAX_FINISHING_SPLITS >= MAX_SPLITS && MAX_FINISHING_SPLITS < (8 * sizeof(unsigned)));
      BOOST_STATIC_ASSERT(log_mean_bin_size >= 0);
      BOOST_STATIC_ASSERT(log_finishing_count >= 0);
      //if we can complete in one iteration, do so
      //This first check allows the compiler to optimize never-executed code out
      if(log_finishing_count < min_size) {
        if(log_range <= min_size && log_range <= MAX_SPLITS) {
          //Return no smaller than a certain minimum limit
          if(log_range <= log_finishing_count)
            return typed_one << log_finishing_count;
          return typed_one << log_range;
        }
      }
      const unsigned base_iterations = MAX_SPLITS - log_min_split_count;
      //sum of n to n + x = ((x + 1) * (n + (n + x)))/2 + log_mean_bin_size
      const unsigned base_range = 
          ((base_iterations + 1) * (MAX_SPLITS + log_min_split_count))/2 
          + log_mean_bin_size;
      //Calculating the required number of iterations, and returning 
      //1 << (iteration_count + min_size)
      if(log_range < base_range) {
        unsigned result = log_min_split_count;
        for(unsigned offset = min_size; offset < log_range; 
          offset += ++result);
        //Preventing overflow; this situation shouldn't occur
        if((result + log_mean_bin_size) >= (8 * sizeof(size_t)))
          return typed_one << ((8 * sizeof(size_t)) - 1);
        return typed_one << (result + log_mean_bin_size);
      }
      //A quick division can calculate the worst-case runtime for larger ranges
      unsigned remainder = log_range - base_range;
      //the MAX_SPLITS - 1 is used to calculate the ceiling of the division
      unsigned bit_length = ((((MAX_SPLITS - 1) + remainder)/MAX_SPLITS) 
        + base_iterations + min_size);
      //Preventing overflow; this situation shouldn't occur
      if(bit_length >= (8 * sizeof(size_t)))
        return typed_one << ((8 * sizeof(size_t)) - 1);
      //n(log_range)/MAX_SPLITS + C, optimizing worst-case performance
      return typed_one << bit_length;
    }

    //Find the minimum and maximum using <
    template <class RandomAccessIter>
    inline void 
    find_extremes(RandomAccessIter current, RandomAccessIter last,
                  RandomAccessIter & max, RandomAccessIter & min)
    {
      min = max = current;
	    //It is assumed we have more than 1 element; there are multiple checks for this
	    while(!(*(current + 1) < *current)) {
        //If everything is in sorted order, return
		    if(++current == last - 1)
			    return;
	    }
	    //The maximum is the last sorted element
	    max = current;
      //Start from the first unsorted element
      while(++current < last) {
        if(*max < *current)
          max = current;
        else if(*current < *min)
          min = current;
      }
    }

    //Uses a user-defined comparison operator to find minimum and maximum
    template <class RandomAccessIter, class Compare>
    inline void 
    find_extremes(RandomAccessIter current, RandomAccessIter last,
                RandomAccessIter & max, RandomAccessIter & min, Compare comp)
    {
      min = max = current;
      while(!comp(*(current + 1), *current)) {
        //If everything is in sorted order, return
		    if(++current == last - 1)
			    return;
	    }
	    //The maximum is the last sorted element
	    max = current;
      while(++current < last) {
        if(comp(*max, *current))
          max = current;
        else if(comp(*current, *min))
          min = current;
      }
    }

    //Gets a non-negative right bit shift to operate as a logarithmic divisor
    template<unsigned log_mean_bin_size>
    inline int
    get_log_divisor(size_t count, unsigned log_range)
    {
      int log_divisor;
      //If we can finish in one iteration without exceeding either 
      //(2 to the MAX_FINISHING_SPLITS) or n bins, do so
      if((log_divisor = log_range - rough_log_2_size(count)) <= 0 
         && log_range <= MAX_FINISHING_SPLITS)
        log_divisor = 0;
      else {
        //otherwise divide the data into an optimized number of pieces
        log_divisor += log_mean_bin_size;
        if(log_divisor < 0)
          log_divisor = 0;
        //Cannot exceed MAX_SPLITS or cache misses slow down bin lookups
        if((log_range - log_divisor) > MAX_SPLITS)
          log_divisor = log_range - MAX_SPLITS;
      }
      return log_divisor;
    }

    template <class RandomAccessIter>
    inline RandomAccessIter * 
    size_bins(std::vector<size_t> &bin_sizes, std::vector<RandomAccessIter>
  &bin_cache, unsigned cache_offset, unsigned &cache_end, unsigned bin_count)
    {
      //Assure space for the size of each bin, followed by initializing sizes
      if(bin_count > bin_sizes.size())
        bin_sizes.resize(bin_count);
      for(size_t u = 0; u < bin_count; u++)
        bin_sizes[u] = 0;
      //Make sure there is space for the bins
      cache_end = cache_offset + bin_count;
      if(cache_end > bin_cache.size())
        bin_cache.resize(cache_end);
      return &(bin_cache[cache_offset]);
    }

    //Implementation for recursive integer sorting
    template <class RandomAccessIter, class Div_type, class Size_type>
    inline void 
    spread_sort_rec(RandomAccessIter first, RandomAccessIter last,
              std::vector<RandomAccessIter> &bin_cache, unsigned cache_offset
              , std::vector<size_t> &bin_sizes)
    {
      //This step is roughly 10% of runtime, but it helps avoid worst-case 
      //behavior and improve behavior with real data
      //If you know the maximum and minimum ahead of time, you can pass those 
      //values in and skip this step for the first iteration
      RandomAccessIter max, min;
      find_extremes(first, last, max, min);
      //max and min will be the same iff all values are equivalent
      if(max == min)
        return;
      RandomAccessIter * target_bin;
      unsigned log_divisor = get_log_divisor<LOG_MEAN_BIN_SIZE>(last - first,
                        rough_log_2_size(Size_type((*max >> 0) - (*min >> 0))));
      Div_type div_min = *min >> log_divisor;
      Div_type div_max = *max >> log_divisor;
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = 
        size_bins(bin_sizes, bin_cache, cache_offset, cache_end, bin_count);
    
      //Calculating the size of each bin; this takes roughly 10% of runtime
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[size_t((*(current++) >> log_divisor) - div_min)]++;
      //Assign the bin positions
      bins[0] = first;
      for(unsigned u = 0; u < bin_count - 1; u++)
        bins[u + 1] = bins[u] + bin_sizes[u];

      RandomAccessIter nextbinstart = first;
      //Swap into place
      //This dominates runtime, mostly in the swap and bin lookups
      for(unsigned u = 0; u < bin_count - 1; ++u) {
        RandomAccessIter * local_bin = bins + u;
        nextbinstart += bin_sizes[u];
        //Iterating over each element in this bin
        for(RandomAccessIter current = *local_bin; current < nextbinstart; 
            ++current) {
          //Swapping elements in current into place until the correct 
          //element has been swapped in
          for(target_bin = (bins + ((*current >> log_divisor) - div_min));  
              target_bin != local_bin; 
            target_bin = bins + ((*current >> log_divisor) - div_min)) {
            //3-way swap; this is about 1% faster than a 2-way swap
            //The main advantage is less copies are involved per item 
            //put in the correct place
            typename std::iterator_traits<RandomAccessIter>::value_type tmp;
            RandomAccessIter b = (*target_bin)++;
            RandomAccessIter * b_bin = bins + ((*b >> log_divisor) - div_min);
            if (b_bin != local_bin) {
              RandomAccessIter c = (*b_bin)++;
              tmp = *c;
              *c = *b;
            } 
            else
              tmp = *b;
            *b = *current;
            *current = tmp;
          }
        }
        *local_bin = nextbinstart;
      }
      bins[bin_count - 1] = last;
  
      //If we've bucketsorted, the array is sorted and we should skip recursion
      if(!log_divisor)
        return;
      //log_divisor is the remaining range; calculating the comparison threshold
      size_t max_count = 
        get_max_count<LOG_MEAN_BIN_SIZE, LOG_MIN_SPLIT_COUNT, LOG_FINISHING_COUNT>(log_divisor);
  
      //Recursing
      RandomAccessIter lastPos = first;
      for(unsigned u = cache_offset; u < cache_end; lastPos = bin_cache[u], 
          ++u) {
        Size_type count = bin_cache[u] - lastPos;
        //don't sort unless there are at least two items to Compare
        if(count < 2)
          continue;
        //using std::sort if its worst-case is better
        if(count < max_count)
          std::sort(lastPos, bin_cache[u]);
        else
          spread_sort_rec<RandomAccessIter, Div_type, Size_type>(lastPos, bin_cache[u],
                                            bin_cache, cache_end, bin_sizes);
      }
    }

    //Generic bitshift-based 3-way swapping code
    template <class RandomAccessIter, class Div_type, class Right_shift>
    inline void inner_swap_loop(RandomAccessIter * bins, 
      const RandomAccessIter & nextbinstart, unsigned ii, Right_shift &rshift
      , const unsigned log_divisor, const Div_type div_min) 
    {
      RandomAccessIter * local_bin = bins + ii;
      for(RandomAccessIter current = *local_bin; current < nextbinstart; 
          ++current) {
        for(RandomAccessIter * target_bin = 
            (bins + (rshift(*current, log_divisor) - div_min));  
            target_bin != local_bin; 
            target_bin = bins + (rshift(*current, log_divisor) - div_min)) {
          typename std::iterator_traits<RandomAccessIter>::value_type tmp;
          RandomAccessIter b = (*target_bin)++;
          RandomAccessIter * b_bin = 
            bins + (rshift(*b, log_divisor) - div_min);
          //Three-way swap; if the item to be swapped doesn't belong 
          //in the current bin, swap it to where it belongs
          if (b_bin != local_bin) {
            RandomAccessIter c = (*b_bin)++;
            tmp = *c;
            *c = *b;
          } 
          //Note: we could increment current once the swap is done in this case
          //but that seems to impair performance
          else
            tmp = *b;
          *b = *current;
          *current = tmp;
        }
      }
      *local_bin = nextbinstart;
    }

    //Standard swapping wrapper for ascending values
    template <class RandomAccessIter, class Div_type, class Right_shift>
    inline void swap_loop(RandomAccessIter * bins, 
             RandomAccessIter & nextbinstart, unsigned ii, Right_shift &rshift
             , const std::vector<size_t> &bin_sizes
             , const unsigned log_divisor, const Div_type div_min) 
    {
      nextbinstart += bin_sizes[ii];
      inner_swap_loop<RandomAccessIter, Div_type, Right_shift>(bins, 
                              nextbinstart, ii, rshift, log_divisor, div_min);
    }

    //Functor implementation for recursive sorting
    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Compare, class Size_type, unsigned log_mean_bin_size, 
                unsigned log_min_split_count, unsigned log_finishing_count>
    inline void 
    spread_sort_rec(RandomAccessIter first, RandomAccessIter last,
          std::vector<RandomAccessIter> &bin_cache, unsigned cache_offset
          , std::vector<size_t> &bin_sizes, Right_shift rshift, Compare comp)
    {
      RandomAccessIter max, min;
      find_extremes(first, last, max, min, comp);
      if(max == min)
        return;
      unsigned log_divisor = get_log_divisor<log_mean_bin_size>(last - first,
            rough_log_2_size(Size_type(rshift(*max, 0) - rshift(*min, 0))));
      Div_type div_min = rshift(*min, log_divisor);
      Div_type div_max = rshift(*max, log_divisor);
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, bin_count);
        
      //Calculating the size of each bin
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[unsigned(rshift(*(current++), log_divisor) - div_min)]++;
      bins[0] = first;
      for(unsigned u = 0; u < bin_count - 1; u++)
        bins[u + 1] = bins[u] + bin_sizes[u];
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      for(unsigned u = 0; u < bin_count - 1; ++u)
        swap_loop<RandomAccessIter, Div_type, Right_shift>(bins, nextbinstart,
                                  u, rshift, bin_sizes, log_divisor, div_min);
      bins[bin_count - 1] = last;
      
      //If we've bucketsorted, the array is sorted
      if(!log_divisor)
        return;
      
      //Recursing
      size_t max_count = get_max_count<log_mean_bin_size, log_min_split_count, 
                          log_finishing_count>(log_divisor);
      RandomAccessIter lastPos = first;
      for(unsigned u = cache_offset; u < cache_end; lastPos = bin_cache[u], 
          ++u) {
        size_t count = bin_cache[u] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[u], comp);
        else
          spread_sort_rec<RandomAccessIter, Div_type, Right_shift, Compare, 
        Size_type, log_mean_bin_size, log_min_split_count, log_finishing_count>
      (lastPos, bin_cache[u], bin_cache, cache_end, bin_sizes, rshift, comp);
      }
    }

    //Functor implementation for recursive sorting with only Shift overridden
    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Size_type, unsigned log_mean_bin_size, 
              unsigned log_min_split_count, unsigned log_finishing_count>
    inline void 
    spread_sort_rec(RandomAccessIter first, RandomAccessIter last,
              std::vector<RandomAccessIter> &bin_cache, unsigned cache_offset
              , std::vector<size_t> &bin_sizes, Right_shift rshift)
    {
      RandomAccessIter max, min;
      find_extremes(first, last, max, min);
      if(max == min)
        return;
      unsigned log_divisor = get_log_divisor<log_mean_bin_size>(last - first,
            rough_log_2_size(Size_type(rshift(*max, 0) - rshift(*min, 0))));
      Div_type div_min = rshift(*min, log_divisor);
      Div_type div_max = rshift(*max, log_divisor);
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, bin_count);
        
      //Calculating the size of each bin
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[unsigned(rshift(*(current++), log_divisor) - div_min)]++;
      bins[0] = first;
      for(unsigned u = 0; u < bin_count - 1; u++)
        bins[u + 1] = bins[u] + bin_sizes[u];
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      for(unsigned ii = 0; ii < bin_count - 1; ++ii)
        swap_loop<RandomAccessIter, Div_type, Right_shift>(bins, nextbinstart,
                                ii, rshift, bin_sizes, log_divisor, div_min);
      bins[bin_count - 1] = last;
      
      //If we've bucketsorted, the array is sorted
      if(!log_divisor)
        return;
      
      //Recursing
      size_t max_count = get_max_count<log_mean_bin_size, log_min_split_count, 
                          log_finishing_count>(log_divisor);
      RandomAccessIter lastPos = first;
      for(unsigned u = cache_offset; u < cache_end; lastPos = bin_cache[u], 
          ++u) {
        size_t count = bin_cache[u] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[u]);
        else
          spread_sort_rec<RandomAccessIter, Div_type, Right_shift, Size_type,
          log_mean_bin_size, log_min_split_count, log_finishing_count>(lastPos,
                      bin_cache[u], bin_cache, cache_end, bin_sizes, rshift);
      }
    }

    //Holds the bin vector and makes the initial recursive call
    template <class RandomAccessIter, class Div_type>
    //Only use spreadsort if the integer can fit in a size_t
    inline typename boost::enable_if_c< sizeof(Div_type) <= sizeof(size_t), void >::type
    integer_sort(RandomAccessIter first, RandomAccessIter last, Div_type)
    {
        std::vector<size_t> bin_sizes;
        std::vector<RandomAccessIter> bin_cache;
        spread_sort_rec<RandomAccessIter, Div_type, size_t>(first, last, 
          bin_cache, 0, bin_sizes);
    }

    //Holds the bin vector and makes the initial recursive call
    template <class RandomAccessIter, class Div_type>
    //Only use spreadsort if the integer can fit in a uintmax_t
    inline typename boost::enable_if_c< (sizeof(Div_type) > sizeof(size_t)) 
      && sizeof(Div_type) <= sizeof(boost::uintmax_t), void >::type
    integer_sort(RandomAccessIter first, RandomAccessIter last, Div_type)
    {
        std::vector<size_t> bin_sizes;
        std::vector<RandomAccessIter> bin_cache;
        spread_sort_rec<RandomAccessIter, Div_type, boost::uintmax_t>(first, 
          last, bin_cache, 0, bin_sizes);
    }

    template <class RandomAccessIter, class Div_type>
    inline typename boost::disable_if_c< sizeof(Div_type) <= sizeof(size_t) 
      || sizeof(Div_type) <= sizeof(boost::uintmax_t), void >::type
    //defaulting to std::sort when integer_sort won't work
    integer_sort(RandomAccessIter first, RandomAccessIter last, Div_type)
    {
      //Warning that we're using std::sort, even though integer_sort was called
      BOOST_STATIC_WARNING( sizeof(Div_type) <= sizeof(size_t) );
      std::sort(first, last);
    }

    
    //Same for the full functor version
    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Compare>
    //Only use spreadsort if the integer can fit in a size_t
    inline typename boost::enable_if_c< sizeof(Div_type) <= sizeof(size_t), 
                                 void >::type
    integer_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
                Right_shift shift, Compare comp)
    {
        std::vector<size_t> bin_sizes;
        std::vector<RandomAccessIter> bin_cache;
        spread_sort_rec<RandomAccessIter, Div_type, Right_shift, Compare, 
          size_t, LOG_MEAN_BIN_SIZE, LOG_MIN_SPLIT_COUNT, LOG_FINISHING_COUNT>
          (first, last, bin_cache, 0, bin_sizes, shift, comp);
    }
   
    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Compare>
    //Only use spreadsort if the integer can fit in a uintmax_t
    inline typename boost::enable_if_c< (sizeof(Div_type) > sizeof(size_t)) 
      && sizeof(Div_type) <= sizeof(boost::uintmax_t), void >::type
    integer_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
                Right_shift shift, Compare comp)
    {
        std::vector<size_t> bin_sizes;
        std::vector<RandomAccessIter> bin_cache;
        spread_sort_rec<RandomAccessIter, Div_type, Right_shift, Compare, 
                        boost::uintmax_t, LOG_MEAN_BIN_SIZE, 
                        LOG_MIN_SPLIT_COUNT, LOG_FINISHING_COUNT>
          (first, last, bin_cache, 0, bin_sizes, shift, comp);
    }

    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Compare>
    inline typename boost::disable_if_c< sizeof(Div_type) <= sizeof(size_t) 
      || sizeof(Div_type) <= sizeof(boost::uintmax_t), void >::type
    //defaulting to std::sort when integer_sort won't work
    integer_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
                Right_shift shift, Compare comp)
    {
      //Warning that we're using std::sort, even though integer_sort was called
      BOOST_STATIC_WARNING( sizeof(Div_type) <= sizeof(size_t) );
      std::sort(first, last, comp);
    }

              
    //Same for the right shift version
    template <class RandomAccessIter, class Div_type, class Right_shift>
    //Only use spreadsort if the integer can fit in a size_t
    inline typename boost::enable_if_c< sizeof(Div_type) <= sizeof(size_t), 
                                 void >::type
    integer_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
                Right_shift shift)
    {
        std::vector<size_t> bin_sizes;
        std::vector<RandomAccessIter> bin_cache;
        spread_sort_rec<RandomAccessIter, Div_type, Right_shift, size_t, 
          LOG_MEAN_BIN_SIZE, LOG_MIN_SPLIT_COUNT, LOG_FINISHING_COUNT>
          (first, last, bin_cache, 0, bin_sizes, shift);
    }

    template <class RandomAccessIter, class Div_type, class Right_shift>
    //Only use spreadsort if the integer can fit in a uintmax_t
    inline typename boost::enable_if_c< (sizeof(Div_type) > sizeof(size_t)) 
      && sizeof(Div_type) <= sizeof(boost::uintmax_t), void >::type
    integer_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
                Right_shift shift)
    {
        std::vector<size_t> bin_sizes;
        std::vector<RandomAccessIter> bin_cache;
        spread_sort_rec<RandomAccessIter, Div_type, Right_shift, 
                        boost::uintmax_t, LOG_MEAN_BIN_SIZE, 
                        LOG_MIN_SPLIT_COUNT, LOG_FINISHING_COUNT>
          (first, last, bin_cache, 0, bin_sizes, shift);
    }

    template <class RandomAccessIter, class Div_type, class Right_shift>
    inline typename boost::disable_if_c< sizeof(Div_type) <= sizeof(size_t) 
      || sizeof(Div_type) <= sizeof(boost::uintmax_t), void >::type
    //defaulting to std::sort when integer_sort won't work
    integer_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
                Right_shift shift)
    {
      //Warning that we're using std::sort, even though integer_sort was called
      BOOST_STATIC_WARNING( sizeof(Div_type) <= sizeof(size_t) );
      std::sort(first, last);
    }

    //------------------------------------------------- float_sort details
 
    //Casts a RandomAccessIter to the specified integer type
    template<class Cast_type, class RandomAccessIter>
    inline Cast_type
    cast_float_iter(const RandomAccessIter & floatiter)
    {
      typedef typename std::iterator_traits<RandomAccessIter>::value_type
        Data_type;
      //Only cast IEEE floating-point numbers, and only to same-sized integers
      BOOST_STATIC_ASSERT(sizeof(Cast_type) == sizeof(Data_type));
      BOOST_STATIC_ASSERT(std::numeric_limits<Data_type>::is_iec559);
      BOOST_STATIC_ASSERT(std::numeric_limits<Cast_type>::is_integer);
      Cast_type result;
      std::memcpy(&result, &(*floatiter), sizeof(Data_type));
      return result;
    }

    template <class RandomAccessIter, class Div_type, class Right_shift>
    inline void 
    find_extremes(RandomAccessIter current, RandomAccessIter last, 
                  Div_type & max, Div_type & min, Right_shift rshift)
    {
      min = max = rshift(*current, 0);
      while(++current < last) {
        Div_type value = rshift(*current, 0);
        if(max < value)
          max = value;
        else if(value < min)
          min = value;
      }
    }

    //Specialized swap loops for floating-point casting
    template <class RandomAccessIter, class Div_type>
    inline void inner_float_swap_loop(RandomAccessIter * bins, 
                        const RandomAccessIter & nextbinstart, unsigned ii
                        , const unsigned log_divisor, const Div_type div_min) 
    {
      RandomAccessIter * local_bin = bins + ii;
      for(RandomAccessIter current = *local_bin; current < nextbinstart; 
          ++current) {
        for(RandomAccessIter * target_bin = 
            (bins + ((cast_float_iter<Div_type, RandomAccessIter>(current) >>
                      log_divisor) - div_min));  target_bin != local_bin; 
          target_bin = bins + ((cast_float_iter<Div_type, RandomAccessIter>
                               (current) >> log_divisor) - div_min)) {
          typename std::iterator_traits<RandomAccessIter>::value_type tmp;
          RandomAccessIter b = (*target_bin)++;
          RandomAccessIter * b_bin = bins + ((cast_float_iter<Div_type,
                              RandomAccessIter>(b) >> log_divisor) - div_min);
          //Three-way swap; if the item to be swapped doesn't belong in the 
          //current bin, swap it to where it belongs
          if (b_bin != local_bin) {
            RandomAccessIter c = (*b_bin)++;
            tmp = *c;
            *c = *b;
          } 
          else
            tmp = *b;
          *b = *current;
          *current = tmp;
        }
      }
      *local_bin = nextbinstart;
    }

    template <class RandomAccessIter, class Div_type>
    inline void float_swap_loop(RandomAccessIter * bins, 
                          RandomAccessIter & nextbinstart, unsigned ii, 
                          const std::vector<size_t> &bin_sizes, 
                          const unsigned log_divisor, const Div_type div_min) 
    {
      nextbinstart += bin_sizes[ii];
      inner_float_swap_loop<RandomAccessIter, Div_type>
        (bins, nextbinstart, ii, log_divisor, div_min);
    }

    template <class RandomAccessIter, class Cast_type>
    inline void 
    find_extremes(RandomAccessIter current, RandomAccessIter last, 
                  Cast_type & max, Cast_type & min)
    {
      min = max = cast_float_iter<Cast_type, RandomAccessIter>(current);
      while(++current < last) {
        Cast_type value = cast_float_iter<Cast_type, RandomAccessIter>(current);
        if(max < value)
          max = value;
        else if(value < min)
          min = value;
      }
    }

    //Special-case sorting of positive floats with casting
    template <class RandomAccessIter, class Div_type, class Size_type>
    inline void 
    positive_float_sort_rec(RandomAccessIter first, RandomAccessIter last,
              std::vector<RandomAccessIter> &bin_cache, unsigned cache_offset
              , std::vector<size_t> &bin_sizes)
    {
      Div_type max, min;
      find_extremes<RandomAccessIter, Div_type>(first, last, max, min);
      if(max == min)
        return;
      unsigned log_divisor = 
        get_log_divisor<FLOAT_LOG_MEAN_BIN_SIZE>(last - first, rough_log_2_size(Size_type(max - min)));
      Div_type div_min = min >> log_divisor;
      Div_type div_max = max >> log_divisor;
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, bin_count);
        
      //Calculating the size of each bin
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[unsigned((cast_float_iter<Div_type, RandomAccessIter>(current++) >>
                   log_divisor) - div_min)]++;
      bins[0] = first;
      for(unsigned u = 0; u < bin_count - 1; u++)
        bins[u + 1] = bins[u] + bin_sizes[u];

      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      for(unsigned u = 0; u < bin_count - 1; ++u)
        float_swap_loop<RandomAccessIter, Div_type>
          (bins, nextbinstart, u, bin_sizes, log_divisor, div_min);
      bins[bin_count - 1] = last;
      
      //Return if we've completed bucketsorting
      if(!log_divisor)
        return;
      
      //Recursing
      size_t max_count = get_max_count<FLOAT_LOG_MEAN_BIN_SIZE, FLOAT_LOG_MIN_SPLIT_COUNT, FLOAT_LOG_FINISHING_COUNT>(log_divisor);
      RandomAccessIter lastPos = first;
      for(unsigned u = cache_offset; u < cache_end; lastPos = bin_cache[u], 
          ++u) {
        size_t count = bin_cache[u] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[u]);
        else
          positive_float_sort_rec<RandomAccessIter, Div_type, Size_type>
            (lastPos, bin_cache[u], bin_cache, cache_end, bin_sizes);
      }
    }

    //Sorting negative floats
    //Bins are iterated in reverse because max_neg_float = min_neg_int
    template <class RandomAccessIter, class Div_type, class Size_type>
    inline void 
    negative_float_sort_rec(RandomAccessIter first, RandomAccessIter last,
                        std::vector<RandomAccessIter> &bin_cache, 
                        unsigned cache_offset, std::vector<size_t> &bin_sizes)
    {
      Div_type max, min;
      find_extremes<RandomAccessIter, Div_type>(first, last, max, min);
      if(max == min)
        return;
      unsigned log_divisor = 
        get_log_divisor<FLOAT_LOG_MEAN_BIN_SIZE>(last - first, rough_log_2_size(Size_type(max - min)));
      Div_type div_min = min >> log_divisor;
      Div_type div_max = max >> log_divisor;
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, bin_count);
        
      //Calculating the size of each bin
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[unsigned((cast_float_iter<Div_type, RandomAccessIter>(current++) >>
                   log_divisor) - div_min)]++;
      bins[bin_count - 1] = first;
      for(int ii = bin_count - 2; ii >= 0; --ii)
        bins[ii] = bins[ii + 1] + bin_sizes[ii + 1];
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      //The last bin will always have the correct elements in it
      for(int ii = bin_count - 1; ii > 0; --ii)
        float_swap_loop<RandomAccessIter, Div_type>
          (bins, nextbinstart, ii, bin_sizes, log_divisor, div_min);
      //Update the end position because we don't process the last bin
      bin_cache[cache_offset] = last;
      
      //Return if we've completed bucketsorting
      if(!log_divisor)
        return;
      
      //Recursing
      size_t max_count = get_max_count<FLOAT_LOG_MEAN_BIN_SIZE, FLOAT_LOG_MIN_SPLIT_COUNT, FLOAT_LOG_FINISHING_COUNT>(log_divisor);
      RandomAccessIter lastPos = first;
      for(int ii = cache_end - 1; ii >= (int)cache_offset; 
          lastPos = bin_cache[ii], --ii) {
        size_t count = bin_cache[ii] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[ii]);
        else
          negative_float_sort_rec<RandomAccessIter, Div_type, Size_type>
            (lastPos, bin_cache[ii], bin_cache, cache_end, bin_sizes);
      }
    }

    //Sorting negative floats
    //Bins are iterated in reverse order because max_neg_float = min_neg_int
    template <class RandomAccessIter, class Div_type, class Right_shift, class Size_type>
    inline void 
    negative_float_sort_rec(RandomAccessIter first, RandomAccessIter last,
              std::vector<RandomAccessIter> &bin_cache, unsigned cache_offset
              , std::vector<size_t> &bin_sizes, Right_shift rshift)
    {
      Div_type max, min;
      find_extremes(first, last, max, min, rshift);
      if(max == min)
        return;
      unsigned log_divisor = 
        get_log_divisor<FLOAT_LOG_MEAN_BIN_SIZE>(last - first, rough_log_2_size(Size_type(max - min)));
      Div_type div_min = min >> log_divisor;
      Div_type div_max = max >> log_divisor;
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, bin_count);
        
      //Calculating the size of each bin
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[unsigned(rshift(*(current++), log_divisor) - div_min)]++;
      bins[bin_count - 1] = first;
      for(int ii = bin_count - 2; ii >= 0; --ii)
        bins[ii] = bins[ii + 1] + bin_sizes[ii + 1];
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      //The last bin will always have the correct elements in it
      for(int ii = bin_count - 1; ii > 0; --ii)
        swap_loop<RandomAccessIter, Div_type, Right_shift>
          (bins, nextbinstart, ii, rshift, bin_sizes, log_divisor, div_min);
      //Update the end position of the unprocessed last bin
      bin_cache[cache_offset] = last;
      
      //Return if we've completed bucketsorting
      if(!log_divisor)
        return;
      
      //Recursing
      size_t max_count = get_max_count<FLOAT_LOG_MEAN_BIN_SIZE, FLOAT_LOG_MIN_SPLIT_COUNT, FLOAT_LOG_FINISHING_COUNT>(log_divisor);
      RandomAccessIter lastPos = first;
      for(int ii = cache_end - 1; ii >= (int)cache_offset; 
          lastPos = bin_cache[ii], --ii) {
        size_t count = bin_cache[ii] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[ii]);
        else
          negative_float_sort_rec<RandomAccessIter, Div_type, Right_shift, Size_type>
            (lastPos, bin_cache[ii], bin_cache, cache_end, bin_sizes, rshift);
      }
    }

    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Compare, class Size_type>
    inline void 
    negative_float_sort_rec(RandomAccessIter first, RandomAccessIter last,
            std::vector<RandomAccessIter> &bin_cache, unsigned cache_offset,
            std::vector<size_t> &bin_sizes, Right_shift rshift, Compare comp)
    {
      Div_type max, min;
      find_extremes(first, last, max, min, rshift);
      if(max == min)
        return;
      unsigned log_divisor = 
        get_log_divisor<FLOAT_LOG_MEAN_BIN_SIZE>(last - first, rough_log_2_size(Size_type(max - min)));
      Div_type div_min = min >> log_divisor;
      Div_type div_max = max >> log_divisor;
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, bin_count);
        
      //Calculating the size of each bin
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[unsigned(rshift(*(current++), log_divisor) - div_min)]++;
      bins[bin_count - 1] = first;
      for(int ii = bin_count - 2; ii >= 0; --ii)
        bins[ii] = bins[ii + 1] + bin_sizes[ii + 1];
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      //The last bin will always have the correct elements in it
      for(int ii = bin_count - 1; ii > 0; --ii)
        swap_loop<RandomAccessIter, Div_type, Right_shift>
          (bins, nextbinstart, ii, rshift, bin_sizes, log_divisor, div_min);
      //Update the end position of the unprocessed last bin
      bin_cache[cache_offset] = last;
      
      //Return if we've completed bucketsorting
      if(!log_divisor)
        return;
      
      //Recursing
      size_t max_count = get_max_count<FLOAT_LOG_MEAN_BIN_SIZE, FLOAT_LOG_MIN_SPLIT_COUNT, FLOAT_LOG_FINISHING_COUNT>(log_divisor);
      RandomAccessIter lastPos = first;
      for(int ii = cache_end - 1; ii >= (int)cache_offset; 
          lastPos = bin_cache[ii], --ii) {
        size_t count = bin_cache[ii] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[ii], comp);
        else
          negative_float_sort_rec<RandomAccessIter, Div_type, Right_shift,
                                  Compare, Size_type>(lastPos, bin_cache[ii], bin_cache,
                                          cache_end, bin_sizes, rshift, comp);
      }
    }

    //Casting special-case for floating-point sorting
    template <class RandomAccessIter, class Div_type, class Size_type>
    inline void 
    float_sort_rec(RandomAccessIter first, RandomAccessIter last,
                std::vector<RandomAccessIter> &bin_cache, unsigned cache_offset
                , std::vector<size_t> &bin_sizes)
    {
      Div_type max, min;
      find_extremes<RandomAccessIter, Div_type>(first, last, max, min);
      if(max == min)
        return;
      unsigned log_divisor = 
        get_log_divisor<FLOAT_LOG_MEAN_BIN_SIZE>(last - first, rough_log_2_size(Size_type(max - min)));
      Div_type div_min = min >> log_divisor;
      Div_type div_max = max >> log_divisor;
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, bin_count);
        
      //Calculating the size of each bin
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[unsigned((cast_float_iter<Div_type, RandomAccessIter>(current++) >>
                   log_divisor) - div_min)]++;
      //The index of the first positive bin
      //Must be divided small enough to fit into an integer
      unsigned first_positive = (div_min < 0) ? unsigned(-div_min) : 0;
      //Resetting if all bins are negative
      if(cache_offset + first_positive > cache_end)
        first_positive = cache_end - cache_offset;
      //Reversing the order of the negative bins
      //Note that because of the negative/positive ordering direction flip
      //We can not depend upon bin order and positions matching up
      //so bin_sizes must be reused to contain the end of the bin
      if(first_positive > 0) {
        bins[first_positive - 1] = first;
        for(int ii = first_positive - 2; ii >= 0; --ii) {
          bins[ii] = first + bin_sizes[ii + 1];
          bin_sizes[ii] += bin_sizes[ii + 1];
        }
        //Handling positives following negatives
        if(first_positive < bin_count) {
          bins[first_positive] = first + bin_sizes[0];
          bin_sizes[first_positive] += bin_sizes[0];
        }
      }
      else
        bins[0] = first;
      for(unsigned u = first_positive; u < bin_count - 1; u++) {
        bins[u + 1] = first + bin_sizes[u];
        bin_sizes[u + 1] += bin_sizes[u];
      }
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      for(unsigned u = 0; u < bin_count; ++u) {
        nextbinstart = first + bin_sizes[u];
        inner_float_swap_loop<RandomAccessIter, Div_type>
          (bins, nextbinstart, u, log_divisor, div_min);
      }
      
      if(!log_divisor)
        return;
      
      //Handling negative values first
      size_t max_count = get_max_count<FLOAT_LOG_MEAN_BIN_SIZE, FLOAT_LOG_MIN_SPLIT_COUNT, FLOAT_LOG_FINISHING_COUNT>(log_divisor);
      RandomAccessIter lastPos = first;
      for(int ii = cache_offset + first_positive - 1; ii >= (int)cache_offset;
          lastPos = bin_cache[ii--]) {
        size_t count = bin_cache[ii] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[ii]);
        //sort negative values using reversed-bin spread_sort
        else
          negative_float_sort_rec<RandomAccessIter, Div_type, Size_type>
            (lastPos, bin_cache[ii], bin_cache, cache_end, bin_sizes);
      }
      
      for(unsigned u = cache_offset + first_positive; u < cache_end; 
          lastPos = bin_cache[u], ++u) {
        size_t count = bin_cache[u] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[u]);
        //sort positive values using normal spread_sort
        else
          positive_float_sort_rec<RandomAccessIter, Div_type, Size_type>
            (lastPos, bin_cache[u], bin_cache, cache_end, bin_sizes);
      }
    }

    //Functor implementation for recursive sorting
    template <class RandomAccessIter, class Div_type, class Right_shift
      , class Size_type>
    inline void 
    float_sort_rec(RandomAccessIter first, RandomAccessIter last,
              std::vector<RandomAccessIter> &bin_cache, unsigned cache_offset
              , std::vector<size_t> &bin_sizes, Right_shift rshift)
    {
      Div_type max, min;
      find_extremes(first, last, max, min, rshift);
      if(max == min)
        return;
      unsigned log_divisor = 
        get_log_divisor<FLOAT_LOG_MEAN_BIN_SIZE>(last - first, rough_log_2_size(Size_type(max - min)));
      Div_type div_min = min >> log_divisor;
      Div_type div_max = max >> log_divisor;
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, bin_count);
        
      //Calculating the size of each bin
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[unsigned(rshift(*(current++), log_divisor) - div_min)]++;
      //The index of the first positive bin
      unsigned first_positive = (div_min < 0) ? unsigned(-div_min) : 0;
      //Resetting if all bins are negative
      if(cache_offset + first_positive > cache_end)
        first_positive = cache_end - cache_offset;
      //Reversing the order of the negative bins
      //Note that because of the negative/positive ordering direction flip
      //We can not depend upon bin order and positions matching up
      //so bin_sizes must be reused to contain the end of the bin
      if(first_positive > 0) {
        bins[first_positive - 1] = first;
        for(int ii = first_positive - 2; ii >= 0; --ii) {
          bins[ii] = first + bin_sizes[ii + 1];
          bin_sizes[ii] += bin_sizes[ii + 1];
        }
        //Handling positives following negatives
        if((unsigned)first_positive < bin_count) {
          bins[first_positive] = first + bin_sizes[0];
          bin_sizes[first_positive] += bin_sizes[0];
        }
      }
      else
        bins[0] = first;
      for(unsigned u = first_positive; u < bin_count - 1; u++) {
        bins[u + 1] = first + bin_sizes[u];
        bin_sizes[u + 1] += bin_sizes[u];
      }
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      for(unsigned u = 0; u < bin_count; ++u) {
        nextbinstart = first + bin_sizes[u];
        inner_swap_loop<RandomAccessIter, Div_type, Right_shift>
          (bins, nextbinstart, u, rshift, log_divisor, div_min);
      }
      
      //Return if we've completed bucketsorting
      if(!log_divisor)
        return;
      
      //Handling negative values first
      size_t max_count = get_max_count<FLOAT_LOG_MEAN_BIN_SIZE, FLOAT_LOG_MIN_SPLIT_COUNT, FLOAT_LOG_FINISHING_COUNT>(log_divisor);
      RandomAccessIter lastPos = first;
      for(int ii = cache_offset + first_positive - 1; ii >= (int)cache_offset;
          lastPos = bin_cache[ii--]) {
        size_t count = bin_cache[ii] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[ii]);
        //sort negative values using reversed-bin spread_sort
        else
          negative_float_sort_rec<RandomAccessIter, Div_type,
            Right_shift, Size_type>(lastPos, bin_cache[ii], bin_cache, cache_end,
                         bin_sizes, rshift);
      }
      
      for(unsigned u = cache_offset + first_positive; u < cache_end; 
          lastPos = bin_cache[u], ++u) {
        size_t count = bin_cache[u] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[u]);
        //sort positive values using normal spread_sort
        else
          spread_sort_rec<RandomAccessIter, Div_type, Right_shift, Size_type, 
                          FLOAT_LOG_MEAN_BIN_SIZE, FLOAT_LOG_MIN_SPLIT_COUNT, 
                          FLOAT_LOG_FINISHING_COUNT>
            (lastPos, bin_cache[u], bin_cache, cache_end, bin_sizes, rshift);
      }
    }

    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Compare, class Size_type>
    inline void 
    float_sort_rec(RandomAccessIter first, RandomAccessIter last,
            std::vector<RandomAccessIter> &bin_cache, unsigned cache_offset,
            std::vector<size_t> &bin_sizes, Right_shift rshift, Compare comp)
    {
      Div_type max, min;
      find_extremes(first, last, max, min, rshift);
      if(max == min)
        return;
      unsigned log_divisor = 
        get_log_divisor<FLOAT_LOG_MEAN_BIN_SIZE>(last - first, rough_log_2_size(Size_type(max - min)));
      Div_type div_min = min >> log_divisor;
      Div_type div_max = max >> log_divisor;
      unsigned bin_count = unsigned(div_max - div_min) + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, bin_count);
        
      //Calculating the size of each bin
      for (RandomAccessIter current = first; current != last;)
        bin_sizes[unsigned(rshift(*(current++), log_divisor) - div_min)]++;
      //The index of the first positive bin
      unsigned first_positive = (div_min < 0) ? unsigned(-div_min) : 0;
      //Resetting if all bins are negative
      if(cache_offset + first_positive > cache_end)
        first_positive = cache_end - cache_offset;
      //Reversing the order of the negative bins
      //Note that because of the negative/positive ordering direction flip
      //We can not depend upon bin order and positions matching up
      //so bin_sizes must be reused to contain the end of the bin
      if(first_positive > 0) {
        bins[first_positive - 1] = first;
        for(int ii = first_positive - 2; ii >= 0; --ii) {
          bins[ii] = first + bin_sizes[ii + 1];
          bin_sizes[ii] += bin_sizes[ii + 1];
        }
        //Handling positives following negatives
        if((unsigned)first_positive < bin_count) {
          bins[first_positive] = first + bin_sizes[0];
          bin_sizes[first_positive] += bin_sizes[0];
        }
      }
      else
        bins[0] = first;
      for(unsigned u = first_positive; u < bin_count - 1; u++) {
        bins[u + 1] = first + bin_sizes[u];
        bin_sizes[u + 1] += bin_sizes[u];
      }
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      for(unsigned u = 0; u < bin_count; ++u) {
        nextbinstart = first + bin_sizes[u];
        inner_swap_loop<RandomAccessIter, Div_type, Right_shift>
          (bins, nextbinstart, u, rshift, log_divisor, div_min);
      }
      
      //Return if we've completed bucketsorting
      if(!log_divisor)
        return;
      
      //Handling negative values first
      size_t max_count = get_max_count<FLOAT_LOG_MEAN_BIN_SIZE, FLOAT_LOG_MIN_SPLIT_COUNT, FLOAT_LOG_FINISHING_COUNT>(log_divisor);
      RandomAccessIter lastPos = first;
      for(int ii = cache_offset + first_positive - 1; ii >= (int)cache_offset;
          lastPos = bin_cache[ii--]) {
        size_t count = bin_cache[ii] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[ii], comp);
        //sort negative values using reversed-bin spread_sort
        else
          negative_float_sort_rec<RandomAccessIter, Div_type,
                                  Right_shift, Compare, Size_type>(lastPos, bin_cache[ii],
                              bin_cache, cache_end, bin_sizes, rshift, comp);
      }
      
      for(unsigned u = cache_offset + first_positive; u < cache_end; 
          lastPos = bin_cache[u], ++u) {
        size_t count = bin_cache[u] - lastPos;
        if(count < 2)
          continue;
        if(count < max_count)
          std::sort(lastPos, bin_cache[u], comp);
        //sort positive values using normal spread_sort
        else
          spread_sort_rec<RandomAccessIter, Div_type, Right_shift, Compare, 
                          Size_type, FLOAT_LOG_MEAN_BIN_SIZE, 
                          FLOAT_LOG_MIN_SPLIT_COUNT, FLOAT_LOG_FINISHING_COUNT>
      (lastPos, bin_cache[u], bin_cache, cache_end, bin_sizes, rshift, comp);
      }
    }

    //Checking whether the value type is a float, and trying a 32-bit integer
    template <class RandomAccessIter>
    inline typename boost::enable_if_c< sizeof(boost::uint32_t) == 
      sizeof(typename std::iterator_traits<RandomAccessIter>::value_type) 
      && std::numeric_limits<typename 
      std::iterator_traits<RandomAccessIter>::value_type>::is_iec559, 
      void >::type
    float_sort(RandomAccessIter first, RandomAccessIter last)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      float_sort_rec<RandomAccessIter, boost::int32_t, boost::uint32_t>
        (first, last, bin_cache, 0, bin_sizes);
    }

    //Checking whether the value type is a double, and using a 64-bit integer
    template <class RandomAccessIter>
    inline typename boost::enable_if_c< sizeof(boost::uint64_t) == 
      sizeof(typename std::iterator_traits<RandomAccessIter>::value_type) 
      && std::numeric_limits<typename 
      std::iterator_traits<RandomAccessIter>::value_type>::is_iec559, 
      void >::type
    float_sort(RandomAccessIter first, RandomAccessIter last)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      float_sort_rec<RandomAccessIter, boost::int64_t, boost::uint64_t>
        (first, last, bin_cache, 0, bin_sizes);
    }

    template <class RandomAccessIter>
    inline typename boost::disable_if_c< (sizeof(boost::uint64_t) == 
      sizeof(typename std::iterator_traits<RandomAccessIter>::value_type)
      || sizeof(boost::uint32_t) == 
      sizeof(typename std::iterator_traits<RandomAccessIter>::value_type))
      && std::numeric_limits<typename 
      std::iterator_traits<RandomAccessIter>::value_type>::is_iec559, 
      void >::type
    float_sort(RandomAccessIter first, RandomAccessIter last)
    {
      BOOST_STATIC_WARNING(!(sizeof(boost::uint64_t) == 
      sizeof(typename std::iterator_traits<RandomAccessIter>::value_type)
      || sizeof(boost::uint32_t) == 
      sizeof(typename std::iterator_traits<RandomAccessIter>::value_type))
      || !std::numeric_limits<typename 
      std::iterator_traits<RandomAccessIter>::value_type>::is_iec559);
      std::sort(first, last);
    }

    //These approaches require the user to do the typecast
    //with rshift but default comparision
    template <class RandomAccessIter, class Div_type, class Right_shift>
    inline typename boost::enable_if_c< sizeof(size_t) >= sizeof(Div_type), 
      void >::type
    float_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
               Right_shift rshift)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      float_sort_rec<RandomAccessIter, Div_type, Right_shift, size_t>
        (first, last, bin_cache, 0, bin_sizes, rshift);
    }

    //maximum integer size with rshift but default comparision
    template <class RandomAccessIter, class Div_type, class Right_shift>
    inline typename boost::enable_if_c< sizeof(size_t) < sizeof(Div_type) 
      && sizeof(boost::uintmax_t) >= sizeof(Div_type), void >::type
    float_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
               Right_shift rshift)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      float_sort_rec<RandomAccessIter, Div_type, Right_shift, boost::uintmax_t>
        (first, last, bin_cache, 0, bin_sizes, rshift);
    }

    //sizeof(Div_type) doesn't match, so use std::sort
    template <class RandomAccessIter, class Div_type, class Right_shift>
    inline typename boost::disable_if_c< sizeof(boost::uintmax_t) >= sizeof(Div_type),
      void >::type
    float_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
               Right_shift rshift)
    {
      BOOST_STATIC_WARNING(sizeof(boost::uintmax_t) >= sizeof(Div_type));
      std::sort(first, last);
    }

    //specialized comparison
    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Compare>
    inline typename boost::enable_if_c< sizeof(size_t) >= sizeof(Div_type), 
      void >::type
    float_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
               Right_shift rshift, Compare comp)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      float_sort_rec<RandomAccessIter, Div_type, Right_shift, Compare,
        size_t>
        (first, last, bin_cache, 0, bin_sizes, rshift, comp);
    }

    //max-sized integer with specialized comparison
    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Compare>
    inline typename boost::enable_if_c< sizeof(size_t) < sizeof(Div_type) 
      && sizeof(boost::uintmax_t) >= sizeof(Div_type), void >::type
    float_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
               Right_shift rshift, Compare comp)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      float_sort_rec<RandomAccessIter, Div_type, Right_shift, Compare, 
        boost::uintmax_t>
        (first, last, bin_cache, 0, bin_sizes, rshift, comp);
    }

    //sizeof(Div_type) doesn't match, so use std::sort
    template <class RandomAccessIter, class Div_type, class Right_shift, 
              class Compare>
    inline typename boost::disable_if_c< sizeof(boost::uintmax_t) >= sizeof(Div_type),
      void >::type
    float_sort(RandomAccessIter first, RandomAccessIter last, Div_type,
               Right_shift rshift, Compare comp)
    {
      BOOST_STATIC_WARNING(sizeof(boost::uintmax_t) >= sizeof(Div_type));
      std::sort(first, last, comp);
    }

  //------------------------------------------------- string_sort details

    //Offsetting on identical characters.  This function works a character 
    //at a time for optimal worst-case performance.
    template<class RandomAccessIter>
    inline void
    update_offset(RandomAccessIter first, RandomAccessIter finish, 
                  unsigned &char_offset)
    {
      unsigned nextOffset = char_offset;
      bool done = false;
      while(!done) {
        RandomAccessIter curr = first;
        do {
          //ignore empties, but if the nextOffset would exceed the length or 
          //not match, exit; we've found the last matching character
          if((*curr).size() > char_offset && ((*curr).size() <= 
           (nextOffset + 1) || (*curr)[nextOffset] != (*first)[nextOffset])) {
            done = true;
            break;
          }
        } while(++curr != finish);
        if(!done)
          ++nextOffset;
      } 
      char_offset = nextOffset;
    }

    //Offsetting on identical characters.  This function works a character 
    //at a time for optimal worst-case performance.
    template<class RandomAccessIter, class Get_char, class Get_length>
    inline void
    update_offset(RandomAccessIter first, RandomAccessIter finish, 
                  unsigned &char_offset, Get_char getchar, Get_length length)
    {
      unsigned nextOffset = char_offset;
      bool done = false;
      while(!done) {
        RandomAccessIter curr = first;
        do {
          //ignore empties, but if the nextOffset would exceed the length or 
          //not match, exit; we've found the last matching character
          if(length(*curr) > char_offset && (length(*curr) <= (nextOffset + 1)
          || getchar((*curr), nextOffset) != getchar((*first), nextOffset))) {
            done = true;
            break;
          }
        } while(++curr != finish);
        if(!done)
          ++nextOffset;
      } 
      char_offset = nextOffset;
    }

    //This comparison functor assumes strings are identical up to char_offset
    template<class Data_type, class Unsigned_char_type>
    struct offset_lessthan {
      offset_lessthan(unsigned char_offset) : fchar_offset(char_offset){}
      inline bool operator()(const Data_type &x, const Data_type &y) const 
      {
        unsigned minSize = std::min(x.size(), y.size());
        for(unsigned u = fchar_offset; u < minSize; ++u) {
          BOOST_STATIC_ASSERT(sizeof(x[u]) == sizeof(Unsigned_char_type));
          if(static_cast<Unsigned_char_type>(x[u]) <
             static_cast<Unsigned_char_type>(y[u]))
            return true;
          else if(static_cast<Unsigned_char_type>(y[u]) <
                  static_cast<Unsigned_char_type>(x[u]))
            return false;
        }
        return x.size() < y.size();
      }
      unsigned fchar_offset;
    };

    //Compares strings assuming they are identical up to char_offset
    template<class Data_type, class Unsigned_char_type>
    struct offset_greaterthan {
      offset_greaterthan(unsigned char_offset) : fchar_offset(char_offset){}
      inline bool operator()(const Data_type &x, const Data_type &y) const 
      {
        unsigned minSize = std::min(x.size(), y.size());
        for(unsigned u = fchar_offset; u < minSize; ++u) {
          BOOST_STATIC_ASSERT(sizeof(x[u]) == sizeof(Unsigned_char_type));
          if(static_cast<Unsigned_char_type>(x[u]) >
             static_cast<Unsigned_char_type>(y[u]))
            return true;
          else if(static_cast<Unsigned_char_type>(y[u]) >
                  static_cast<Unsigned_char_type>(x[u]))
            return false;
        }
        return x.size() > y.size();
      }
      unsigned fchar_offset;
    };

    //This comparison functor assumes strings are identical up to char_offset
    template<class Data_type, class Get_char, class Get_length>
    struct offset_char_lessthan {
      offset_char_lessthan(unsigned char_offset) : fchar_offset(char_offset){}
      inline bool operator()(const Data_type &x, const Data_type &y) const 
      {
        unsigned minSize = std::min(length(x), length(y));
        for(unsigned u = fchar_offset; u < minSize; ++u) {
          if(getchar(x, u) < getchar(y, u))
            return true;
          else if(getchar(y, u) < getchar(x, u))
            return false;
        }
        return length(x) < length(y);
      }
      unsigned fchar_offset;
      Get_char getchar;
      Get_length length;
    };

    //String sorting recursive implementation
    template <class RandomAccessIter, class Unsigned_char_type>
    inline void 
    string_sort_rec(RandomAccessIter first, RandomAccessIter last, 
                unsigned char_offset, std::vector<RandomAccessIter> &bin_cache
                , unsigned cache_offset, std::vector<size_t> &bin_sizes)
    {
      typedef typename std::iterator_traits<RandomAccessIter>::value_type
        Data_type;
      //This section makes handling of long identical substrings much faster
      //with a mild average performance impact.
      //Iterate to the end of the empties.  If all empty, return
      while((*first).size() <= char_offset) {
        if(++first == last)
          return;
      }
      RandomAccessIter finish = last - 1;
      //Getting the last non-empty
      for(;(*finish).size() <= char_offset; --finish);
      ++finish;
      //Offsetting on identical characters.  This section works 
      //a character at a time for optimal worst-case performance.
      update_offset(first, finish, char_offset);
      
      const unsigned bin_count = (1 << (sizeof(Unsigned_char_type)*8));
      //Equal worst-case of radix and comparison is when bin_count = n*log(n).
      const unsigned max_size = bin_count;
      const unsigned membin_count = bin_count + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, membin_count) + 1;
        
      //Calculating the size of each bin; this takes roughly 10% of runtime
      for (RandomAccessIter current = first; current != last; ++current) {
        if((*current).size() <= char_offset) {
          bin_sizes[0]++;
        }
        else
          bin_sizes[static_cast<Unsigned_char_type>((*current)[char_offset]) 
                    + 1]++;
      }
      //Assign the bin positions
      bin_cache[cache_offset] = first;
      for(unsigned u = 0; u < membin_count - 1; u++)
        bin_cache[cache_offset + u + 1] = 
          bin_cache[cache_offset + u] + bin_sizes[u];
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      //handling empty bins
      RandomAccessIter * local_bin = &(bin_cache[cache_offset]);
      nextbinstart +=  bin_sizes[0];
      RandomAccessIter * target_bin;
      //Iterating over each element in the bin of empties
      for(RandomAccessIter current = *local_bin; current < nextbinstart; 
          ++current) {
        //empties belong in this bin
        while((*current).size() > char_offset) {
          target_bin = 
            bins + static_cast<Unsigned_char_type>((*current)[char_offset]);
          iter_swap(current, (*target_bin)++);
        }
      }
      *local_bin = nextbinstart;
      //iterate backwards to find the last bin with elements in it
      //this saves iterations in multiple loops
      unsigned last_bin = bin_count - 1;
      for(; last_bin && !bin_sizes[last_bin + 1]; --last_bin);
      //This dominates runtime, mostly in the swap and bin lookups
      for(unsigned u = 0; u < last_bin; ++u) {
        local_bin = bins + u;
        nextbinstart += bin_sizes[u + 1];
        //Iterating over each element in this bin
        for(RandomAccessIter current = *local_bin; current < nextbinstart; 
            ++current) {
          //Swapping into place until the correct element has been swapped in
          for(target_bin = bins + static_cast<Unsigned_char_type>
              ((*current)[char_offset]);  target_bin != local_bin;
            target_bin = bins + static_cast<Unsigned_char_type>
              ((*current)[char_offset])) iter_swap(current, (*target_bin)++);
        }
        *local_bin = nextbinstart;
      }
      bins[last_bin] = last;
      //Recursing
      RandomAccessIter lastPos = bin_cache[cache_offset];
      //Skip this loop for empties
      for(unsigned u = cache_offset + 1; u < cache_offset + last_bin + 2; 
          lastPos = bin_cache[u], ++u) {
        size_t count = bin_cache[u] - lastPos;
        //don't sort unless there are at least two items to Compare
        if(count < 2)
          continue;
        //using std::sort if its worst-case is better
        if(count < max_size)
          std::sort(lastPos, bin_cache[u], 
              offset_lessthan<Data_type, Unsigned_char_type>(char_offset + 1));
        else
          string_sort_rec<RandomAccessIter, Unsigned_char_type>(lastPos,
              bin_cache[u], char_offset + 1, bin_cache, cache_end, bin_sizes);
      }
    }

    //Sorts strings in reverse order, with empties at the end
    template <class RandomAccessIter, class Unsigned_char_type>
    inline void 
    reverse_string_sort_rec(RandomAccessIter first, RandomAccessIter last,
                unsigned char_offset, std::vector<RandomAccessIter> &bin_cache
                , unsigned cache_offset, std::vector<size_t> &bin_sizes)
    {
      typedef typename std::iterator_traits<RandomAccessIter>::value_type
        Data_type;
      //This section makes handling of long identical substrings much faster
      //with a mild average performance impact.
      RandomAccessIter curr = first;
      //Iterate to the end of the empties.  If all empty, return
      while((*curr).size() <= char_offset) {
        if(++curr == last)
          return;
      }
      //Getting the last non-empty
      while((*(--last)).size() <= char_offset);
      ++last;
      //Offsetting on identical characters.
      update_offset(curr, last, char_offset);
      RandomAccessIter * target_bin;
      
      const unsigned bin_count = (1 << (sizeof(Unsigned_char_type)*8));
      //Equal worst-case of radix and comparison when bin_count = n*log(n).
      const unsigned max_size = bin_count;
      const unsigned membin_count = bin_count + 1;
      const unsigned max_bin = bin_count - 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, membin_count);
      RandomAccessIter * end_bin = &(bin_cache[cache_offset + max_bin]);
        
      //Calculating the size of each bin; this takes roughly 10% of runtime
      for (RandomAccessIter current = first; current != last; ++current) {
        if((*current).size() <= char_offset) {
          bin_sizes[bin_count]++;
        }
        else
          bin_sizes[max_bin - static_cast<Unsigned_char_type>
            ((*current)[char_offset])]++;
      }
      //Assign the bin positions
      bin_cache[cache_offset] = first;
      for(unsigned u = 0; u < membin_count - 1; u++)
        bin_cache[cache_offset + u + 1] = 
          bin_cache[cache_offset + u] + bin_sizes[u];
      
      //Swap into place
      RandomAccessIter nextbinstart = last;
      //handling empty bins
      RandomAccessIter * local_bin = &(bin_cache[cache_offset + bin_count]);
      RandomAccessIter lastFull = *local_bin;
      //Iterating over each element in the bin of empties
      for(RandomAccessIter current = *local_bin; current < nextbinstart; 
          ++current) {
        //empties belong in this bin
        while((*current).size() > char_offset) {
          target_bin = 
            end_bin - static_cast<Unsigned_char_type>((*current)[char_offset]);
          iter_swap(current, (*target_bin)++);
        }
      }
      *local_bin = nextbinstart;
      nextbinstart = first;
      //iterate backwards to find the last non-empty bin 
      //this saves iterations in multiple loops
      unsigned last_bin = max_bin;
      for(; last_bin && !bin_sizes[last_bin]; --last_bin);
      //This dominates runtime, mostly in the swap and bin lookups
      for(unsigned u = 0; u < last_bin; ++u) {
        local_bin = bins + u;
        nextbinstart += bin_sizes[u];
        //Iterating over each element in this bin
        for(RandomAccessIter current = *local_bin; current < nextbinstart; 
            ++current) {
          //Swapping into place until the correct element has been swapped in
          for(target_bin = 
            end_bin - static_cast<Unsigned_char_type>((*current)[char_offset]); 
            target_bin != local_bin; 
            target_bin = 
            end_bin - static_cast<Unsigned_char_type>((*current)[char_offset]))
              iter_swap(current, (*target_bin)++);
        }
        *local_bin = nextbinstart;
      }
      bins[last_bin] = lastFull;
      //Recursing
      RandomAccessIter lastPos = first;
      //Skip this loop for empties
      for(unsigned u = cache_offset; u <= cache_offset + last_bin; 
          lastPos = bin_cache[u], ++u) {
        size_t count = bin_cache[u] - lastPos;
        //don't sort unless there are at least two items to Compare
        if(count < 2)
          continue;
        //using std::sort if its worst-case is better
        if(count < max_size)
          std::sort(lastPos, bin_cache[u], offset_greaterthan<Data_type,
                    Unsigned_char_type>(char_offset + 1));
        else
          reverse_string_sort_rec<RandomAccessIter, Unsigned_char_type>
    (lastPos, bin_cache[u], char_offset + 1, bin_cache, cache_end, bin_sizes);
      }
    }

    //String sorting recursive implementation
    template <class RandomAccessIter, class Unsigned_char_type, class Get_char,
              class Get_length>
    inline void 
    string_sort_rec(RandomAccessIter first, RandomAccessIter last, 
              unsigned char_offset, std::vector<RandomAccessIter> &bin_cache,
              unsigned cache_offset, std::vector<size_t> &bin_sizes, 
              Get_char getchar, Get_length length)
    {
      typedef typename std::iterator_traits<RandomAccessIter>::value_type 
        Data_type;
      //This section makes handling of long identical substrings much faster
      //with a mild average performance impact.
      //Iterate to the end of the empties.  If all empty, return
      while(length(*first) <= char_offset) {
        if(++first == last)
          return;
      }
      RandomAccessIter finish = last - 1;
      //Getting the last non-empty
      for(;length(*finish) <= char_offset; --finish);
      ++finish;
      update_offset(first, finish, char_offset, getchar, length);
      
      const unsigned bin_count = (1 << (sizeof(Unsigned_char_type)*8));
      //Equal worst-case of radix and comparison is when bin_count = n*log(n).
      const unsigned max_size = bin_count;
      const unsigned membin_count = bin_count + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, membin_count) + 1;
        
      //Calculating the size of each bin; this takes roughly 10% of runtime
      for (RandomAccessIter current = first; current != last; ++current) {
        if(length(*current) <= char_offset) {
          bin_sizes[0]++;
        }
        else
          bin_sizes[getchar((*current), char_offset) + 1]++;
      }
      //Assign the bin positions
      bin_cache[cache_offset] = first;
      for(unsigned u = 0; u < membin_count - 1; u++)
        bin_cache[cache_offset + u + 1] = 
          bin_cache[cache_offset + u] + bin_sizes[u];
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      //handling empty bins
      RandomAccessIter * local_bin = &(bin_cache[cache_offset]);
      nextbinstart +=  bin_sizes[0];
      RandomAccessIter * target_bin;
      //Iterating over each element in the bin of empties
      for(RandomAccessIter current = *local_bin; current < nextbinstart; 
          ++current) {
        //empties belong in this bin
        while(length(*current) > char_offset) {
          target_bin = bins + getchar((*current), char_offset);
          iter_swap(current, (*target_bin)++);
        }
      }
      *local_bin = nextbinstart;
      //iterate backwards to find the last bin with elements in it 
      //this saves iterations in multiple loops
      unsigned last_bin = bin_count - 1;
      for(; last_bin && !bin_sizes[last_bin + 1]; --last_bin);
      //This dominates runtime, mostly in the swap and bin lookups
      for(unsigned ii = 0; ii < last_bin; ++ii) {
        local_bin = bins + ii;
        nextbinstart += bin_sizes[ii + 1];
        //Iterating over each element in this bin
        for(RandomAccessIter current = *local_bin; current < nextbinstart; 
            ++current) {
          //Swapping into place until the correct element has been swapped in
          for(target_bin = bins + getchar((*current), char_offset);  
              target_bin != local_bin; 
              target_bin = bins + getchar((*current), char_offset))
            iter_swap(current, (*target_bin)++);
        }
        *local_bin = nextbinstart;
      }
      bins[last_bin] = last;
      
      //Recursing
      RandomAccessIter lastPos = bin_cache[cache_offset];
      //Skip this loop for empties
      for(unsigned u = cache_offset + 1; u < cache_offset + last_bin + 2; 
          lastPos = bin_cache[u], ++u) {
        size_t count = bin_cache[u] - lastPos;
        //don't sort unless there are at least two items to Compare
        if(count < 2)
          continue;
        //using std::sort if its worst-case is better
        if(count < max_size)
          std::sort(lastPos, bin_cache[u], offset_char_lessthan<Data_type,
                    Get_char, Get_length>(char_offset + 1));
        else
          string_sort_rec<RandomAccessIter, Unsigned_char_type, Get_char,
            Get_length>(lastPos, bin_cache[u], char_offset + 1, bin_cache,
                        cache_end, bin_sizes, getchar, length);
      }
    }

    //String sorting recursive implementation
    template <class RandomAccessIter, class Unsigned_char_type, class Get_char,
              class Get_length, class Compare>
    inline void 
    string_sort_rec(RandomAccessIter first, RandomAccessIter last, 
              unsigned char_offset, std::vector<RandomAccessIter> &bin_cache,
              unsigned cache_offset, std::vector<size_t> &bin_sizes, 
              Get_char getchar, Get_length length, Compare comp)
    {
      //This section makes handling of long identical substrings much faster
      //with a mild average performance impact.
      //Iterate to the end of the empties.  If all empty, return
      while(length(*first) <= char_offset) {
        if(++first == last)
          return;
      }
      RandomAccessIter finish = last - 1;
      //Getting the last non-empty
      for(;length(*finish) <= char_offset; --finish);
      ++finish;
      update_offset(first, finish, char_offset, getchar, length);
      
      const unsigned bin_count = (1 << (sizeof(Unsigned_char_type)*8));
      //Equal worst-case of radix and comparison is when bin_count = n*log(n).
      const unsigned max_size = bin_count;
      const unsigned membin_count = bin_count + 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, membin_count) + 1;
        
      //Calculating the size of each bin; this takes roughly 10% of runtime
      for (RandomAccessIter current = first; current != last; ++current) {
        if(length(*current) <= char_offset) {
          bin_sizes[0]++;
        }
        else
          bin_sizes[getchar((*current), char_offset) + 1]++;
      }
      //Assign the bin positions
      bin_cache[cache_offset] = first;
      for(unsigned u = 0; u < membin_count - 1; u++)
        bin_cache[cache_offset + u + 1] = 
          bin_cache[cache_offset + u] + bin_sizes[u];
      
      //Swap into place
      RandomAccessIter nextbinstart = first;
      //handling empty bins
      RandomAccessIter * local_bin = &(bin_cache[cache_offset]);
      nextbinstart +=  bin_sizes[0];
      RandomAccessIter * target_bin;
      //Iterating over each element in the bin of empties
      for(RandomAccessIter current = *local_bin; current < nextbinstart; 
          ++current) {
        //empties belong in this bin
        while(length(*current) > char_offset) {
          target_bin = bins + getchar((*current), char_offset);
          iter_swap(current, (*target_bin)++);
        }
      }
      *local_bin = nextbinstart;
      //iterate backwards to find the last bin with elements in it
      //this saves iterations in multiple loops
      unsigned last_bin = bin_count - 1;
      for(; last_bin && !bin_sizes[last_bin + 1]; --last_bin);
      //This dominates runtime, mostly in the swap and bin lookups
      for(unsigned u = 0; u < last_bin; ++u) {
        local_bin = bins + u;
        nextbinstart += bin_sizes[u + 1];
        //Iterating over each element in this bin
        for(RandomAccessIter current = *local_bin; current < nextbinstart; 
            ++current) {
          //Swapping into place until the correct element has been swapped in
          for(target_bin = bins + getchar((*current), char_offset);  
              target_bin != local_bin; 
              target_bin = bins + getchar((*current), char_offset))
            iter_swap(current, (*target_bin)++);
        }
        *local_bin = nextbinstart;
      }
      bins[last_bin] = last;
      
      //Recursing
      RandomAccessIter lastPos = bin_cache[cache_offset];
      //Skip this loop for empties
      for(unsigned u = cache_offset + 1; u < cache_offset + last_bin + 2; 
          lastPos = bin_cache[u], ++u) {
        size_t count = bin_cache[u] - lastPos;
        //don't sort unless there are at least two items to Compare
        if(count < 2)
          continue;
        //using std::sort if its worst-case is better
        if(count < max_size)
          std::sort(lastPos, bin_cache[u], comp);
        else
          string_sort_rec<RandomAccessIter, Unsigned_char_type, Get_char, 
                          Get_length, Compare>
            (lastPos, bin_cache[u], char_offset + 1, bin_cache, cache_end,
             bin_sizes, getchar, length, comp);
      }
    }

    //Sorts strings in reverse order, with empties at the end
    template <class RandomAccessIter, class Unsigned_char_type, class Get_char,
              class Get_length, class Compare>
    inline void 
    reverse_string_sort_rec(RandomAccessIter first, RandomAccessIter last,
              unsigned char_offset, std::vector<RandomAccessIter> &bin_cache,
              unsigned cache_offset, std::vector<size_t> &bin_sizes, 
              Get_char getchar, Get_length length, Compare comp)
    {
      //This section makes handling of long identical substrings much faster
      //with a mild average performance impact.
      RandomAccessIter curr = first;
      //Iterate to the end of the empties.  If all empty, return
      while(length(*curr) <= char_offset) {
        if(++curr == last)
          return;
      }
      //Getting the last non-empty
      while(length(*(--last)) <= char_offset);
      ++last;
      //Offsetting on identical characters.  This section works 
      //a character at a time for optimal worst-case performance.
      update_offset(curr, last, char_offset, getchar, length);
      
      const unsigned bin_count = (1 << (sizeof(Unsigned_char_type)*8));
      //Equal worst-case of radix and comparison is when bin_count = n*log(n).
      const unsigned max_size = bin_count;
      const unsigned membin_count = bin_count + 1;
      const unsigned max_bin = bin_count - 1;
      unsigned cache_end;
      RandomAccessIter * bins = size_bins(bin_sizes, bin_cache, cache_offset,
                                          cache_end, membin_count);
      RandomAccessIter *end_bin = &(bin_cache[cache_offset + max_bin]);
        
      //Calculating the size of each bin; this takes roughly 10% of runtime
      for (RandomAccessIter current = first; current != last; ++current) {
        if(length(*current) <= char_offset) {
          bin_sizes[bin_count]++;
        }
        else
          bin_sizes[max_bin - getchar((*current), char_offset)]++;
      }
      //Assign the bin positions
      bin_cache[cache_offset] = first;
      for(unsigned u = 0; u < membin_count - 1; u++)
        bin_cache[cache_offset + u + 1] = 
          bin_cache[cache_offset + u] + bin_sizes[u];
      
      //Swap into place
      RandomAccessIter nextbinstart = last;
      //handling empty bins
      RandomAccessIter * local_bin = &(bin_cache[cache_offset + bin_count]);
      RandomAccessIter lastFull = *local_bin;
      RandomAccessIter * target_bin;
      //Iterating over each element in the bin of empties
      for(RandomAccessIter current = *local_bin; current < nextbinstart; 
          ++current) {
        //empties belong in this bin
        while(length(*current) > char_offset) {
          target_bin = end_bin - getchar((*current), char_offset);
          iter_swap(current, (*target_bin)++);
        }
      }
      *local_bin = nextbinstart;
      nextbinstart = first;
      //iterate backwards to find the last bin with elements in it
      //this saves iterations in multiple loops
      unsigned last_bin = max_bin;
      for(; last_bin && !bin_sizes[last_bin]; --last_bin);
      //This dominates runtime, mostly in the swap and bin lookups
      for(unsigned u = 0; u < last_bin; ++u) {
        local_bin = bins + u;
        nextbinstart += bin_sizes[u];
        //Iterating over each element in this bin
        for(RandomAccessIter current = *local_bin; current < nextbinstart; 
            ++current) {
          //Swapping into place until the correct element has been swapped in
          for(target_bin = end_bin - getchar((*current), char_offset);  
              target_bin != local_bin; 
              target_bin = end_bin - getchar((*current), char_offset))
            iter_swap(current, (*target_bin)++);
        }
        *local_bin = nextbinstart;
      }
      bins[last_bin] = lastFull;
      //Recursing
      RandomAccessIter lastPos = first;
      //Skip this loop for empties
      for(unsigned u = cache_offset; u <= cache_offset + last_bin; 
          lastPos = bin_cache[u], ++u) {
        size_t count = bin_cache[u] - lastPos;
        //don't sort unless there are at least two items to Compare
        if(count < 2)
          continue;
        //using std::sort if its worst-case is better
        if(count < max_size)
          std::sort(lastPos, bin_cache[u], comp);
        else
          reverse_string_sort_rec<RandomAccessIter, Unsigned_char_type, 
                                  Get_char, Get_length, Compare>
            (lastPos, bin_cache[u], char_offset + 1, bin_cache, cache_end,
             bin_sizes, getchar, length, comp);
      }
    }

    //Holds the bin vector and makes the initial recursive call
    template <class RandomAccessIter, class Unsigned_char_type>
    inline typename boost::enable_if_c< sizeof(Unsigned_char_type) <= 2, void >::type
    string_sort(RandomAccessIter first, RandomAccessIter last,
                Unsigned_char_type)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      string_sort_rec<RandomAccessIter, Unsigned_char_type>
        (first, last, 0, bin_cache, 0, bin_sizes);
    }

    template <class RandomAccessIter, class Unsigned_char_type>
    inline typename boost::disable_if_c< sizeof(Unsigned_char_type) <= 2, void >::type
    string_sort(RandomAccessIter first, RandomAccessIter last,
                Unsigned_char_type)
    {
      //Warning that we're using std::sort, even though string_sort was called
      BOOST_STATIC_WARNING( sizeof(Unsigned_char_type) <= 2 );
      std::sort(first, last);
    }

    //Holds the bin vector and makes the initial recursive call
    template <class RandomAccessIter, class Unsigned_char_type>
    inline typename boost::enable_if_c< sizeof(Unsigned_char_type) <= 2, void >::type
    reverse_string_sort(RandomAccessIter first, RandomAccessIter last,
                        Unsigned_char_type)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      reverse_string_sort_rec<RandomAccessIter, Unsigned_char_type>
        (first, last, 0, bin_cache, 0, bin_sizes);
    }

    template <class RandomAccessIter, class Unsigned_char_type>
    inline typename boost::disable_if_c< sizeof(Unsigned_char_type) <= 2, void >::type
    reverse_string_sort(RandomAccessIter first, RandomAccessIter last,
                Unsigned_char_type)
    {
      typedef typename std::iterator_traits<RandomAccessIter>::value_type Data_type;
      //Warning that we're using std::sort, even though string_sort was called
      BOOST_STATIC_WARNING( sizeof(Unsigned_char_type) <= 2 );
      std::sort(first, last, std::greater<Data_type>());
    }

    //Holds the bin vector and makes the initial recursive call
    template <class RandomAccessIter, class Get_char, class Get_length, 
              class Unsigned_char_type>
    inline typename boost::enable_if_c< sizeof(Unsigned_char_type) <= 2, void >::type 
    string_sort(RandomAccessIter first, RandomAccessIter last, 
                Get_char getchar, Get_length length, Unsigned_char_type)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      string_sort_rec<RandomAccessIter, Unsigned_char_type, Get_char,
        Get_length>(first, last, 0, bin_cache, 0, bin_sizes, getchar, length);
    }

    template <class RandomAccessIter, class Get_char, class Get_length, 
              class Unsigned_char_type>
    inline typename boost::disable_if_c< sizeof(Unsigned_char_type) <= 2, void >::type 
    string_sort(RandomAccessIter first, RandomAccessIter last, 
                Get_char getchar, Get_length length, Unsigned_char_type)
    {
      //Warning that we're using std::sort, even though string_sort was called
      BOOST_STATIC_WARNING( sizeof(Unsigned_char_type) <= 2 );
      std::sort(first, last);
    }

    //Holds the bin vector and makes the initial recursive call
    template <class RandomAccessIter, class Get_char, class Get_length, 
              class Compare, class Unsigned_char_type>
    inline typename boost::enable_if_c< sizeof(Unsigned_char_type) <= 2, void >::type  
    string_sort(RandomAccessIter first, RandomAccessIter last, 
        Get_char getchar, Get_length length, Compare comp, Unsigned_char_type)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      string_sort_rec<RandomAccessIter, Unsigned_char_type, Get_char
        , Get_length, Compare>
        (first, last, 0, bin_cache, 0, bin_sizes, getchar, length, comp);
    }

    //disable_if_c was refusing to compile, so rewrote to use enable_if_c
    template <class RandomAccessIter, class Get_char, class Get_length, 
              class Compare, class Unsigned_char_type>
    inline typename boost::enable_if_c< (sizeof(Unsigned_char_type) > 2), void >::type  
    string_sort(RandomAccessIter first, RandomAccessIter last, 
        Get_char getchar, Get_length length, Compare comp, Unsigned_char_type)
    {
      //Warning that we're using std::sort, even though string_sort was called
      BOOST_STATIC_WARNING( sizeof(Unsigned_char_type) <= 2 );
      std::sort(first, last, comp);
    }

    //Holds the bin vector and makes the initial recursive call
    template <class RandomAccessIter, class Get_char, class Get_length, 
              class Compare, class Unsigned_char_type>
    inline typename boost::enable_if_c< sizeof(Unsigned_char_type) <= 2, void >::type  
    reverse_string_sort(RandomAccessIter first, RandomAccessIter last, 
        Get_char getchar, Get_length length, Compare comp, Unsigned_char_type)
    {
      std::vector<size_t> bin_sizes;
      std::vector<RandomAccessIter> bin_cache;
      reverse_string_sort_rec<RandomAccessIter, Unsigned_char_type, Get_char,
                              Get_length, Compare>
        (first, last, 0, bin_cache, 0, bin_sizes, getchar, length, comp);
    }

    template <class RandomAccessIter, class Get_char, class Get_length, 
              class Compare, class Unsigned_char_type>
    inline typename boost::disable_if_c< sizeof(Unsigned_char_type) <= 2, void >::type  
    reverse_string_sort(RandomAccessIter first, RandomAccessIter last, 
        Get_char getchar, Get_length length, Compare comp, Unsigned_char_type)
    {
      //Warning that we're using std::sort, even though string_sort was called
      BOOST_STATIC_WARNING( sizeof(Unsigned_char_type) <= 2 );
      std::sort(first, last, comp);
    }
  }
}

#endif
