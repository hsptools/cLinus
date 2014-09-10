#include<stdio.h>
//#include "Atom3d.h"
#include "linus.h"
void reset_conf(LinusProtein *prot)
{
        int i;
        for(i=0;i<prot->num_atoms;i++)
        {
          prot->atoms[i].fp_distance = prot->orig_fp_distance[i];
        }

        ztox(prot->atoms,0,prot->num_atoms);
}

/*
        o *reset_conformation* - restore the conformation of the
        molecule to that specified in the input PDB file from which
        the molecule was first created


        """Restore the conformation of the molecule as existed in
        the PDB file from which the molecule was initially read.
        This function takes no arguments
        """
*/
