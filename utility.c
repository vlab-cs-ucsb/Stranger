//
//  utility.c
//  strangerlib
//
//  Created by Muath Alkhalaf on 9/16/13.
//  Copyright (c) 2013 Muath Alkhalaf. All rights reserved.
//

#include <stdio.h>
#include "utility.h"
#include "stranger.h"
#include "stranger_lib_internal.h"
#include <locale.h>



int intcmpfunc (const void * a, const void * b)
{
    const unsigned *x = a, *y = b;
    if(*x > *y)
        return 1;
    else
        return(*x < *y) ? -1: 0;
}



unsigned roundToNextPow2(unsigned v){
    // compute the next highest power of 2 of 32-bit v
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}


// searches for a state in an ordered list of states using binary search
// inclusive indices
//   0 <= imin when using truncate toward zero divide
//     imid = (imin+imax)/2;
//   imin unrestricted when using truncate toward minus infinity divide
//     imid = (imin+imax)>>1; or
//     imid = (int)floor((imin+imax)/2.0);
bool findStateBS(const unsigned states[], unsigned keyState, unsigned imin, unsigned imax)
{
    // continually narrow search until just one element remains
    while (imin < imax)
    {
        
        int imid = (imax + imin)/2;
        // code must guarantee the interval is reduced at each iteration
        assert(imid < imax);
        // note: 0 <= imin < imax implies imid will always be less than imax
        
        // reduce the search
        if (states[imid] < keyState)
            imin = imid + 1;
        else
            imax = imid;
    }
    // At exit of while:
    //   if A[] is empty, then imax < imin
    //   otherwise imax == imin
    
    // deferred test for equality
    if ((imax == imin) && (states[imin] == keyState))
        return true;
    else
        return false;
}


char *commaprint(unsigned long long n)
{
	static int comma = '\0';
//	static char retbuf[(4*(sizeof(unsigned long * CHAR_BIT + 2)/3/3+1)];
    static char retbuf[30];
	char *p = &retbuf[sizeof(retbuf)-1];
	int i = 0;
    
	if(comma == '\0') {
		struct lconv *lcp = localeconv();
		if(lcp != NULL) {
			if(lcp->thousands_sep != NULL &&
               *lcp->thousands_sep != '\0')
				comma = *lcp->thousands_sep;
			else	comma = ',';
		}
	}
    
	*p = '\0';
    
	do {
		if(i%3 == 0 && i != 0)
			*--p = comma;
		*--p = '0' + n % 10;
		n /= 10;
		i++;
	} while(n != 0);
    
	return p;
}
