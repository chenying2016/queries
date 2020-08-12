//Templated Spreadsort generic hybrid sorting

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

#ifndef BOOST_SPREAD_SORT_H
#define BOOST_SPREAD_SORT_H
#include <algorithm>
#include <vector>
#include <cstring>
#include <limits>
#include <boost/static_assert.hpp>
#include <boost/algorithm/sorting/constants.hpp>
#include <boost/algorithm/sorting/integer_sort.hpp>
#include <boost/algorithm/sorting/float_sort.hpp>

namespace boost {
  //Generic spread_sort call to integer_sort
  template <class RandomAccessIter>
  inline typename boost::enable_if_c< std::numeric_limits< 
    typename std::iterator_traits<RandomAccessIter>::value_type >::is_integer, 
    void >::type
  spread_sort(RandomAccessIter first, RandomAccessIter last) 
  {
    //Don't sort if it's too small to optimize
    if(last - first < detail::MIN_SORT_SIZE)
      std::sort(first, last);
    else
      detail::integer_sort(first, last, *first >> 0);
  }

  //Generic spread_sort call to float_sort
  template <class RandomAccessIter>
  inline typename boost::enable_if_c< !std::numeric_limits< 
    typename std::iterator_traits<RandomAccessIter>::value_type >::is_integer
    && std::numeric_limits< 
    typename std::iterator_traits<RandomAccessIter>::value_type >::is_iec559, 
    void >::type
  spread_sort(RandomAccessIter first, RandomAccessIter last) 
  {
    float_sort(first, last);
  }

  //Fallback to std::sort
  template <class RandomAccessIter>
  inline typename boost::enable_if_c< !std::numeric_limits< 
    typename std::iterator_traits<RandomAccessIter>::value_type >::is_integer&&
    !std::numeric_limits<
    typename std::iterator_traits<RandomAccessIter>::value_type >::is_iec559,
    void >::type
  spread_sort(RandomAccessIter first, RandomAccessIter last)
  {
    std::sort(first, last);
  }
}

#endif
