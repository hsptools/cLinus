#include <stdio.h>
//#include "Atom3d.h"
#include "linus.h"


#include <string.h>

int get_atom_res_id(Atom3d *atoms, int num_atoms, int resid, char *aname)
{
   int i;

    for (i = 0; i < num_atoms; i++)
    {
        if (!strcmp(atoms[i].name,aname) && atoms[i].resnum == resid)
        {
            return i;
        }
   }
   return -1;
}
