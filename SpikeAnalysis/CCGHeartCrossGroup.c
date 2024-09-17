/* CCGHeartCrossGroup.c /*
/* This is a bare-bones C program whos purpose is to compute
   Multi-unit cross correlograms quickly.  Not intended for
   use on its own.  It is designed to be wrapped by a MATLAB
   function.

   Usage =    [CCG,COUNT] = (CCGHeartCrossGroup(TIMES,MARKS,CONDITIONS,BINSIZE,HALFBINS);
   TIMES ( is the name of a binary file containing N doubles giving the spike times
   MARKS is the name of a binary file containing N unsigned ints giving the spike markers.  Don't use zero!
   CONDITIONS binary file containing N unsigned ints giving the condition markers. Will not compute across conditions  Don't use zero!
   BINSIZE is a number giving the size of the ccg bins in TIMES units
   HALFBINS is the number of bins to compute - so there are a total of nBins = 1+2*HALFBINS bins
   
   These should be: double, uint32,uint32,double,uint32

   NB The spikes MUST be sorted.

   CCG contains unsigned ints containing the counts in each bin
   It is like a 3d array, indexed by [nBins*nMarks*Mark1 + nBins*Mark2 + Bin]
   
   PAIRS contains the spike numbers of every pair of spikes included in the ccg.

  
 */


#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define CHECK
#define STRLEN 10000
#define PAIRBLOCKSIZE 1000000000

unsigned int *Pairs;
unsigned int PairCnt, PairSz;
void AddPair(unsigned int n1, unsigned int n2) {
	unsigned int *pui;

	if (PairSz==0) {
/*		mexPrintf("Allocating pair memory\n"); */
		Pairs = mxMalloc(PAIRBLOCKSIZE*sizeof(unsigned int));
		PairSz = PAIRBLOCKSIZE;
		if (!Pairs) mexErrMsgTxt("Could not allocate memory for pairs");
	}
	/* check if array is full, if so add more memory*/
	if(PairCnt>=PairSz) {
/*		mexPrintf("Reallocating pair memory ... ");
		PairSz += PAIRBLOCKSIZE;
		pui = mxRealloc(Pairs, PairSz);
		mexPrintf("got %x\n", pui);
		if (!pui) {
			mxFree(Pairs);

			mexErrMsgTxt("Could not reallocate memory for pairs");
		}
		Pairs = pui;
*/
       mexPrintf("\n Number of pairs %d\n",PairCnt);
		mexErrMsgTxt("Too many pairs");

	}
	Pairs[PairCnt++] = n1;
	Pairs[PairCnt++] = n2;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

	unsigned int nSpikes, nMarks, HalfBins, nBins, i, CountArraySize, CountIndex, nGroups,nCountIndex,nCountArraySize;
	double *Times;
	double BinSize, FurthestEdge;
	unsigned int *Marks, Mark, *Count,*Groups, Group,*NCount;
	unsigned int CenterSpike, Mark1, Mark2, Bin, Group1,Group2,GroupID;
	int SecondSpike; /* we want to let it go negative so we can stop it there... */
	int BinCounter;
	double Time1, Time2;
	char errstr[STRLEN];

	/* global variables are not initialized on each call to mex fn!! */
	PairCnt = 0; PairSz = 0;

	if (nrhs!=5) {
		mexErrMsgTxt("Must have 5 arguments\nBut listen:  You don't want to use this program.\nUse the MATLAB wrapper function CCG instead.\n");
	}

    
	/* get arguments */
	Times = mxGetPr(prhs[0]);
    nSpikes = mxGetNumberOfElements(prhs[0]);
    
    
    if (mxGetNumberOfElements(prhs[1])!=nSpikes) mexErrMsgTxt("Number of marks ~= number of spikes");
    
    
	Marks = (unsigned int *) mxGetPr(prhs[1]);	
	Groups = (unsigned int *) mxGetPr(prhs[2]);
    
	BinSize = mxGetScalar(prhs[3]);
	HalfBins = (unsigned int) mxGetScalar(prhs[4]);
    
	/* derive other constants */
	nBins        = 1+2*HalfBins;
	FurthestEdge = (BinSize * (HalfBins + 0.5));
	
	/* count nMarks */
	nMarks = 0;
    nGroups = 0;
	for(i=0; i<nSpikes; i++) {
		Mark = Marks[i];
        Group = Groups[i];
		if (Mark>nMarks) {
            nMarks = Mark;
        }
      
        if (Group>nGroups){
            nGroups = Group;
        }
        
		if (Mark==0 || Group==0) {
			mexErrMsgTxt("CCGEngine: No zeros allowed in Marks or Groups");
			abort();
		}
	}

    
	/* allocate output array */
	CountArraySize = nMarks * nMarks * nBins * nGroups*nGroups;
     nCountArraySize = nMarks * nGroups * nGroups;
    
	plhs[0] = mxCreateNumericMatrix(CountArraySize, 1, mxUINT32_CLASS, mxREAL);
	Count = (unsigned int *) mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateNumericMatrix(nCountArraySize, 1, mxUINT32_CLASS, mxREAL);
	NCount = (unsigned int *) mxGetPr(plhs[1]);
    

	if (!Times || !Marks || !Count) {
		mexErrMsgTxt("CCGEngine could not allocate memory!\n");
	}

	/* Now the main program .... */


	for(CenterSpike=0; CenterSpike<nSpikes; CenterSpike++) {
		Mark1 = Marks[CenterSpike];
		Time1 = Times[CenterSpike];
        Group1 = Groups[CenterSpike];
        
          
   
                
		 
                for (Group2 = 1;Group2<=nGroups;Group2++){
                 GroupID = (Group1 - 1) * nGroups + Group2 ;
                 nCountIndex = nMarks * (GroupID-1) + Mark1-1 ;
                 NCount[nCountIndex]++;
                }
        
		/* Go back from CenterSpike */
		for(SecondSpike=CenterSpike-1; SecondSpike>=0; SecondSpike--) {
			Time2 = Times[SecondSpike];

			/* check if we have left the interesting region */
			if((Time1 - Time2) > FurthestEdge) break;

			/* calculate bin */
			Group2 = Groups[SecondSpike];
			
          
            Mark2 = Marks[SecondSpike];
            
    
         
            GroupID = (Group1 - 1) * nGroups + Group2 ;
			 Bin = HalfBins + (int)(floor((0.5+(Time2-Time1)/BinSize)));
			
			/* Bin = HalfBins + (int)(rintf((Time2-Time1)/BinSize)); */

			
			CountIndex = nBins*nMarks*nMarks*(GroupID - 1) + nBins*nMarks*(Mark1-1) + nBins*(Mark2-1) + Bin;
#ifdef CHECK
			if (CountIndex<0 || CountIndex >= CountArraySize) {
			
				sprintf(errstr, "ert a: t1 %f t2 %f m1 %d m2 %d Bin %d, index %d out of bounds",
					Time1, Time2, Mark1, Mark2, BinCounter, CountIndex);
				mexErrMsgTxt(errstr);
			}
#endif

			/* increment count */
			Count[CountIndex]++;
		
            }
		

		/* Now do the same thing going forward... */
		for(SecondSpike=CenterSpike+1; SecondSpike<nSpikes; SecondSpike++) {
			Time2 = Times[SecondSpike];

			/* check if we have left the interesting region */
			if ((Time1 - Time2) <= -FurthestEdge) break;

            Group2 = Groups[SecondSpike];
           
            GroupID = (Group1 - 1) * nGroups + Group2 ;
             Mark2 = Marks[SecondSpike];
            
			/* calculate bin */
			Bin = HalfBins + (int)(((0.5+(Time2-Time1)/BinSize)));

			
			CountIndex =  nBins*nMarks*nMarks*(GroupID - 1) + nBins*nMarks*(Mark1-1) + nBins*(Mark2-1) + Bin;

#ifdef CHECK
			if (CountIndex<0 || CountIndex >= CountArraySize) {
				sprintf(errstr, "err b: t1 %f t2 %f m1 %d m2 %d Bin %d, index %d out %d of %d bounds",
					Time1, Time2, Mark1, Mark2, Bin, CountIndex,(Time1 - Time2),FurthestEdge);
				mexErrMsgTxt(errstr);
			}
#endif

			/* increment count */
			Count[CountIndex]++;
			
            

		}
	}




	/* sayonara */
}
