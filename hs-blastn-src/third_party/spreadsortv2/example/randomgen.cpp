// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/ for updates, documentation, and revision history.

#include <stdio.h>
#include "stdlib.h"
#include <fstream>
#include <iostream>

int main(int argc, const char ** argv) {
	//Always seed with the same value, to get the same results
	srand(1);
	//defaults
	unsigned high_shift = 16;
	unsigned low_shift = 16;
	unsigned count = 1000000;
	//Reading in user arguments
	if(argc > 1)
		high_shift = atoi(argv[1]);
	if(argc > 2)
		low_shift = atoi(argv[2]);
	if(argc > 3)
		count = atoi(argv[3]);
	if(high_shift > 16)
		high_shift = 16;
	if(low_shift > 16)
		low_shift = 16;
	std::ofstream ofile;
	ofile.open("input.txt", std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
	if(ofile.bad()) {
		printf("could not open input.txt for writing!\n");
		return 1;
	}
	//buffering file output for speed
	unsigned uDivideFactor = 1000;
	//Skipping buffering for small files
	if(count < uDivideFactor * 100)
		uDivideFactor = count;
	int * pNumbers = (int *) malloc(uDivideFactor * sizeof(int));
	//Generating semirandom numbers
	for(unsigned u = 0; u < count/uDivideFactor; ++u) {
		unsigned i = 0;
		for(; i< uDivideFactor; ++i) {
      //Generating a 32-bit random number
      pNumbers[i] = (rand() % (1 << low_shift)) | ((rand() % (1 << high_shift)) << 16);
      if(16 == low_shift && rand() % 2)
        pNumbers[i] |= 1 << 15;
      //Adding the sign bit
      if(16 == high_shift && rand() % 2)
        pNumbers[i] *= -1;
		}
		ofile.write( (char *) pNumbers, uDivideFactor * 4 );
	}
  ofile.close();
	return 0;
}
