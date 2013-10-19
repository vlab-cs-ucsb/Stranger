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
#include <string.h>


typedef struct UIntArrayList_ {
    unsigned *list;
    size_t size;
    size_t index;
    bool sorted;
} UIntArrayList, *PUIntArrayList;
PUIntArrayList createUIntArrayList(size_t size);
void insertIntoUIntArrayList(PUIntArrayList pUIntArrayList, unsigned elem);
void insertIntoUIntSortedArrayList(PUIntArrayList pUIntArrayList, unsigned elem);
bool deleteFromUIntArrayList(PUIntArrayList pUIntArrayList, unsigned elem);
bool searchUIntArrayList(PUIntArrayList pUIntArrayList, unsigned elem, size_t *pIndex);
bool searchUIntArrayListBS(PUIntArrayList pUIntArrayList, unsigned elem, size_t *pIndex);
void sortUIntArrayList(PUIntArrayList pUIntArrayList);
void freeUIntArrayList(PUIntArrayList pUIntArrayList);


typedef struct StatePair_ {
    unsigned first;
    unsigned second;
    char *escapedChars;
} StatePair, *PStatePair;

typedef struct StatePairArrayList_ {
    PStatePair *list;
    size_t size;
    size_t index;
    bool sorted;
    size_t numOfEscapedChars;
} StatePairArrayList, *PStatePairArrayList;

PStatePairArrayList createStatePairArrayList(size_t size, size_t numOfEscapedChars);
void insertIntoStatePairArrayList(PStatePairArrayList pStatePairArrayList, unsigned first, unsigned second, char escapeChar);
void insertIntoStatePairSortedArrayList(PStatePairArrayList pStatePairArrayList, unsigned first, unsigned second, char escapeChar);
void addEscapeCharToStatePairArrayList(PStatePairArrayList pStatePairArrayList, unsigned first, unsigned second, char escapeChar);
bool deleteFromStatePairArrayList(PStatePairArrayList pStatePairArrayList, unsigned  first, unsigned second);
bool searchStatePairArrayList(PStatePairArrayList pStatePairArrayList, unsigned first, unsigned second, size_t *pIndex);
bool searchStatePairArrayListBS(PStatePairArrayList pStatePairArrayList, unsigned first, unsigned second, size_t *pIndex);
void sortStatePairArrayList(PStatePairArrayList pStatePairArrayList);
void printStatePairArrayList(PStatePairArrayList pStatePairArrayList);
void freeStatePairArrayList(PStatePairArrayList pStatePairArrayList);

unsigned roundToNextPow2(unsigned v);

int intcmpfunc (const void * a, const void * b);

int statePaircmpfunc (const void * a, const void * b);

bool findStateBS(const unsigned states[], unsigned keyState, unsigned imin, unsigned imax);

bool findStatePairsBS(const PStatePair statePairs[], PStatePair keyStatePair, unsigned imin, unsigned imax);

char *commaprint(unsigned long long n);

#endif
