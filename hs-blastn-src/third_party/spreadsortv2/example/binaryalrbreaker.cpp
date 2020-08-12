// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/ for updates, documentation, and revision history.

#include <boost/algorithm/sorting/spread_sort.hpp>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace boost;
using namespace std;

#define DATA_TYPE boost::uint64_t

#define ALR_THRESHOLD 20

const unsigned max_count = ALR_THRESHOLD - 1;
const unsigned bit_shift = detail::rough_log_2_size(max_count) - detail::LOG_MEAN_BIN_SIZE;
const unsigned radix_threshold = detail::rough_log_2_size(max_count) + 1;

const DATA_TYPE typed_one = 1;

void
fill_vector(vector<DATA_TYPE> & input, const DATA_TYPE base_value, unsigned remaining_bits
            , const vector<unsigned> & indices, int index)
{
  if(index < 0) {
    for(unsigned u = 0; u < max_count; ++u)
      input.push_back((base_value << remaining_bits) + (rand() % (1 << remaining_bits)));
  }
  else {
    unsigned shift = indices[index];
    fill_vector(input, (base_value << shift) + ((1 << shift) - 1), 
                remaining_bits - shift, indices, index - 1);
    fill_vector(input, base_value << shift, remaining_bits - shift, indices, 
                index - 1);
  }
}

//Pass in an argument to test std::sort
int main(int argc, const char ** argv) {
  vector<DATA_TYPE> input;
  vector<unsigned> offsets;
  unsigned total_length = sizeof(uint64_t) * 8;
  unsigned bit_length = total_length;
  unsigned bit_offset = bit_shift;
  bit_length -= radix_threshold;
  for(; bit_length >= ++bit_offset; bit_length -= bit_offset)
    offsets.push_back(bit_offset);
  for(int ii = (1 << bit_length) - 1; ii >= 0; --ii)
    fill_vector(input, ii, total_length - bit_length, offsets, offsets.size() - 1);
	//Run multiple loops, if requested
	for(unsigned u = 0; u < 2; ++u) {
    vector<DATA_TYPE> array = input;
		clock_t start, end;
		double elapsed;
		start = clock();
		if(u)
			std::sort(array.begin(), array.end());
		else
			spread_sort(array.begin(), array.end());
		end = clock();
		elapsed = ((double) (end - start));
    if(u)
      printf("std::sort elapsed time %f\n", elapsed / CLOCKS_PER_SEC);
    else
      printf("spreadsort elapsed time %f\n", elapsed / CLOCKS_PER_SEC);
		array.clear();
	}
	return 0;
}
