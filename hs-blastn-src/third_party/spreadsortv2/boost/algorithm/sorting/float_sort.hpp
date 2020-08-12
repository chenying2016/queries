//Templated Spreadsort-based implementation of float_sort and float_mem_cast

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

#ifndef BOOST_FLOAT_SORT_H
#define BOOST_FLOAT_SORT_H
#include <algorithm>
#include <vector>
#include <cstring>
#include <limits>
#include <boost/static_assert.hpp>
#include <boost/algorithm/sorting/constants.hpp>
#include <boost/algorithm/sorting/detail/spread_sort.hpp>

namespace boost {
  //Casts a float to the specified integer type
  template<class Data_type, class Cast_type>
  inline Cast_type
  float_mem_cast(const Data_type & data)
  {
    //Only cast IEEE floating-point numbers, and only to a same-sized integer
    BOOST_STATIC_ASSERT(sizeof(Cast_type) == sizeof(Data_type));
    BOOST_STATIC_ASSERT(std::numeric_limits<Data_type>::is_iec559);
    BOOST_STATIC_ASSERT(std::numeric_limits<Cast_type>::is_integer);
    Cast_type result;
    std::memcpy(&result, &data, sizeof(Cast_type));
    return result;
  }

  //float_sort with casting to the appropriate size
  template <class RandomAccessIter>
  inline void float_sort(RandomAccessIter first, RandomAccessIter last) 
  {
    if(last - first < detail::MIN_SORT_SIZE)
      std::sort(first, last);
    else
      detail::float_sort(first, last);
  }

  //float_sort with functors
  template <class RandomAccessIter, class Right_shift>
  inline void float_sort(RandomAccessIter first, RandomAccessIter last,
                         Right_shift rshift) 
  {
    if(last - first < detail::MIN_SORT_SIZE)
      std::sort(first, last);
    else
      detail::float_sort(first, last, rshift(*first, 0), rshift);
  }

  template <class RandomAccessIter, class Right_shift, class Compare>
  inline void float_sort(RandomAccessIter first, RandomAccessIter last,
                         Right_shift rshift, Compare comp) 
  {
    if(last - first < detail::MIN_SORT_SIZE)
      std::sort(first, last, comp);
    else
      detail::float_sort(first, last, rshift(*first, 0), rshift, comp);
  }
}

#endif
