#include <stdio.h>
#include <stdlib.h>
//#include "Atom3d.h"
#include "linus.h"


#include <string.h>

int get_atom_ind_with_name(Atom3d *atoms, int start, int end, char *aname)
{
    int i;
    
   for (i = start; i < end; i++)
   {
        
        if (!strcmp(atoms[i].name,aname))
        {
            //printf("Start %d End %d i %d Aname %s Atoms name %s\n", start, end, i, aname, atoms[i].name);    
            //printf("Found atom %d!\n",i);
            return i;
        }
   }
   return -1;
}
