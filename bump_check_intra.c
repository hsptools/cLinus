#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Atom3d.h"
//#include "linus.h"



int bump_check_intra(Atom3d *atoms, int *fai, int res_id )
{
    int i1, i2, start, i;

    i1 = fai[res_id];
    i2 = fai[res_id+1];
    start = i1 + 7;

    // if residue is gly or ala no intra bump
    if (start >= i2)
    {
        return 0;
    }

    // check sidechain with carbonyl oxygen
    for (i=start; i<i2; i++)
    {
        return close(&atoms[i1+3], &atoms[i], -1.0e0); //Oatm
    }

    // check sidechain with N-H
    // (for PRO i1+5 is HA and this won't clash any intra)

    for (i=start; i<i2; i++)
    {
        return close(&atoms[i1+5], &atoms[i], -1.0e0); //Hatm
    }

    return 0;
}

