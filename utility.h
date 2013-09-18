//
//  utility.h
//  strangerlib
//
//  Created by Muath Alkhalaf on 9/16/13.
//  Copyright (c) 2013 Muath Alkhalaf. All rights reserved.
//

#ifndef strangerlib_utility_h
#define strangerlib_utility_h

#include <stdbool.h>

unsigned roundToNextPow2(unsigned v);

int intcmpfunc (const void * a, const void * b);

bool findStateBS(const unsigned states[], unsigned keyState, unsigned imin, unsigned imax);

char *commaprint(unsigned long long n);

#endif
