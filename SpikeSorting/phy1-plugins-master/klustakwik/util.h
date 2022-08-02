 /*
 * util.h
 *
 *  Created on: 11 Nov 2011
 *      Author: dan
 */

#ifndef UTIL_H_
#define UTIL_H_

#include<stdio.h>
#include "numerics.h"
#include "log.h"

extern const scalar HugeScore;

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

integer irand(integer min, integer max);
FILE *fopen_safe(char *fname, char *mode);
void MatPrint(FILE *fp, scalar *Mat, integer nRows, integer nCols);
void CompareVectors(scalar *A, scalar *B, integer N);

#endif /* UTIL_H_ */
