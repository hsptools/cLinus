/*
	Combined Multiple Recursive Random Number Generator.

	This implementation is borrowed entirely from the C code in

	L'Ecuyer, Pierre., "Good parameters and implementations for combined
	multiple recursive random number generators,"

	(http://www.iro.umontreal.ca/~lecuyer/myftp/papers/combmrg2.ps)

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NORM 2.328306549295728e-10
#define M1   4294967087.0
#define M2   4294944443.0
#define A12     1403580.0
#define A13N     810728.0
#define A21      527612.0
#define A23N    1370589.0

/*
generator  calls cmrg32_cmrg32
  cmrg32_cmrg32 - Checks if all numbers are double and calls new_cmrg32Object
      new_cmrg32Object - check if all values != 0; 0 < 1,2,3 <M1;
 

*/

double random1(double *s)
{

    /* check that not all seeds of each set are zero */
    
    if(s[0]==0.0 && s[1]==0.0 && s[2]==0.0 && s[3]==0.0 && s[4]==0.0 && s[5]==0.0)
    {
        printf("All seeds are 0\n");
        exit(1);
    }
    
    /* check that the first set of seeds are in the range 0 to M1 */
    
    if ((s[0]<0.0e0 || s[0]>=M1) || (s[1]<0.0e0 || s[1]>=M1) || (s[2]<0.0e0 || s[2]>=M1)) 
    {
        printf("First three seeds must be in range 1 to 4294967087\n");
		exit(1);
	}
    
    /* check that the first set of seeds are in the range 0 to M2 */

	if ((s[3]<0.0e0 || s[3]>=M2) || (s[4]<0.0e0 || s[4]>=M2) || (s[5]<0.0e0 || s[5]>=M2)) 
    {
		printf("Last three seeds must be in range 1 to 4294944443\n");
		exit(1);
	}

    /* we need to ensure that the seeds are exactly integral */
    
    s[0] = floor(s[0]);
    s[1] = floor(s[1]);
    s[2] = floor(s[2]);
    s[3] = floor(s[3]);
    s[4] = floor(s[4]);
    s[5] = floor(s[5]);
    
    long k;
    double p1, p2;
	
    p1 = A12 * s[1] - A13N * s[0];
	k = p1 / M1;
	p1 -= k * M1;
	if (p1 < 0.0)
		p1 += M1;
	s[0] = s[1];
	s[1] = s[2];
	s[2] = p1;

	p2 = A21 * s[5] - A23N * s[3];
	k  = p2 / M2;
	p2 -= k * M2;
	if (p2 < 0.0)
		p2 += M2;
	s[3] = s[4];
	s[4] = s[5];
	s[5] = p2;

    double ran;
    
    
	if (p1 <= p2)
		ran = ((p1 - p2 + M1) * NORM);
	else
		ran = ((p1 - p2) * NORM);
    //printf("Random %lf\n",ran);
    return ran;
}







