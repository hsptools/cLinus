#include <stdio.h>
#include "LinusRes.h"
//#include "Atom3d.h"
#include "linus.h"


int get_scacceptor(Atom3d * atoms, int start, int end, char * resname, int * scatoms)
{
    int s, e; // get residues range from linusRes

    get_res_range(resname,&s,&e); // get indexes for linusRes for a given residue
    //printf("Linusres Res %s start %d end %d\n",resname,s,e);
    int i,ctr=0;

    char *a= " O  ";


    for(i=s;i<e;i++)
    {
       //printf("%s-%s-%s-%d\n",resname,atm_name_arr[i],a,hba_arr[i]);
       if(strcmp(atm_name_arr[i],a) && hba_arr[i])
       {
            scatoms[ctr] = get_atom_ind_with_name(atoms,start,end, atm_name_arr[i]);
            //printf("Acceptor %s %s\n",atoms[scatoms[ctr]].name,atoms[atoms[scatoms[ctr]].fp_ind].name);
            ctr++;
       }
    }
    return ctr;
}

