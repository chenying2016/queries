//Templated hybrid string_sort

//          Copyright Steven J. Ross 2001 - 2009.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org/ for updates, documentation, and revision history
// See http://www.boost.org/libs/algorithm/sorting for library home page.
      
/*
Some improvements suggested by:
Phil Endecott and Frank Gennari
*/

#ifndef BOOST_STRING_SORT_H
#define BOOST_STRING_SORT_H
#include <algorithm>
#include <vector>
#include <cstring>
#include <limits>
#include <boost/static_assert.hpp>
#include <boost/algorithm/sorting/constants.hpp>
#include <boost/algorithm/sorting/detail/spread_sort.hpp>

namespace boost {
  //Allows character-type overloads
  template <class RandomAccessIter, class Unsigned_char_type>
  inline void string_sort(RandomAccessIter first, RandomAccessIter last,
                          Unsigned_char_type unused) 
  {
    //Don't sort if it's too small to optimize
    if(last - first < detail::MIN_SORT_SIZE)
      std::sort(first, last);
    else
      detail::string_sort(first, last, unused);
  }

  //Top-level sorting call; wraps using default of unsigned char
  template <class RandomAccessIter>
  inline void string_sort(RandomAccessIter first, RandomAccessIter last) 
  {
    unsigned char unused = '\0';
    string_sort(first, last, unused);
  }

  //Allows character-type overloads
  template <class RandomAccessIter, class Compare, class Unsigned_char_type>
  inline void reverse_string_sort(RandomAccessIter first, 
                RandomAccessIter last, Compare comp, Unsigned_char_type unused) 
  {
    //Don't sort if it's too small to optimize
    if(last - first < detail::MIN_SORT_SIZE)
      std::sort(first, last, comp);
    else
      detail::reverse_string_sort(first, last, unused);
  }

  //Top-level sorting call; wraps using default of unsigned char
  template <class RandomAccessIter, class Compare>
  inline void reverse_string_sort(RandomAccessIter first, 
                                  RandomAccessIter last, Compare comp) 
  {
    unsigned char unused = '\0';
    reverse_string_sort(first, last, comp, unused);
  }

  template <class RandomAccessIter, class Get_char, class Get_length>
  inline void string_sort(RandomAccessIter first, RandomAccessIter last,
                          Get_char getchar, Get_length length) 
  {
    //Don't sort if it's too small to optimize
    if(last - first < detail::MIN_SORT_SIZE)
      std::sort(first, last);
    else {
      //skipping past empties, which allows us to get the character type 
      //.empty() is not used so as not to require a user declaration of it
      while(!length(*first)) {
        if(++first == last)
          return;
      }
      detail::string_sort(first, last, getchar, length, getchar((*first), 0));
    }
  }

  template <class RandomAccessIter, class Get_char, class Get_length, 
            class Compare>
  inline void string_sort(RandomAccessIter first, RandomAccessIter last,
                          Get_char getchar, Get_length length, Compare comp) 
  {
    //Don't sort if it's too small to optimize
    if(last - first < detail::MIN_SORT_SIZE)
      std::sort(first, last, comp);
    else {
      //skipping past empties, which allows us to get the character type 
      //.empty() is not used so as not to require a user declaration of it
      while(!length(*first)) {
        if(++first == last)
          return;
      }
      detail::string_sort(first, last, getchar, length, comp, 
                          getchar((*first), 0));
    }
  }

  template <class RandomAccessIter, class Get_char, class Get_length, 
            class Compare>
  inline void reverse_string_sort(RandomAccessIter first, 
    RandomAccessIter last, Get_char getchar, Get_length length, Compare comp) 
  {
    //Don't sort if it's too small to optimize
    if(last - first < detail::MIN_SORT_SIZE)
      std::sort(first, last, comp);
    else {
      //skipping past empties, which allows us to get the character type 
      //.empty() is not used so as not to require a user declaration of it
      while(!length(*(--last))) {
        //If there is just one non-empty at the beginning, this is sorted
        if(first == last)
          return;
      }
      //making last just after the end of the non-empty part of the array
      detail::reverse_string_sort(first, last + 1, getchar, length, comp,
                                  getchar((*last), 0));
    }
  }
}

#endif
