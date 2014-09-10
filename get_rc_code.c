#include <stdio.h>
#include <string.h>
#include <math.h>
#include "linusSecStr.h"

void get_rc_code(int phi, int psi, char *code)
{
    printf("Searching rc\n");
    int i;
    for(i=0;i<144;i++)
    {
        if(rc_phi[i] == phi && rc_psi[i] == psi)
        {
            
            strcpy(code,rc_code[i]);
            //printf("getrc index %d Code is %s Saved as %s\n",i,rc_code[i],code);
        }
    }

}

