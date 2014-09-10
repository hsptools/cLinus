#include <stdio.h>
//#include "Atom3d.h"
#include "linus.h"


 /*Construct parameters for accessible surface area scoring

    Arguments

        o *mol* - instance of linusMol - molecule whose accessible surface
          area is to be calculated

        o *hphob* - boolean - whether to include hydrophobic atoms in
          the evaluation of surface area (by default, all atoms are used
          in steric calculations except ACE, NME)

        o *hphil* - boolean - whether to include hydrophillic atoms in
          the evaluation of surface area (by default, all atoms are used
          in steric calculations except ACE, NME)

        o *hydrogens* boolean - whether to include hydrogens in the
          scoring function (as contributing to the return value or
          as steric hinderances to solvation).  If included, hydrogens
          are hydrophobic, according to *HYDROPHILLIC* in **linusRes**.

    Result

        Returns an object that can be passed to *asa_score* or
        *print_asa_score* functions in the **linusScore** module.
    */

int make_asa_list(LinusProtein *p, int hphob, int hphil, int hydrogens, int *flags)
{
    int PRESENT = 1;
    int CONTRIB = 2;
  /*
    // This should be in  Main/Calling function
    int flags[p->num_atoms];

    */
    if(hphob < 0)
        hphob = 1;
    if(hphil < 0)
        hphil = 0;
    if(hydrogens <0)
        hydrogens=0;

    int i;
    for (i=0;i<p->num_atoms;i++)
    {
        flags[i] = 0;
    }
    int minres,maxres;
    get_res_extents(p,&minres,&maxres);

    int r;
    for(r=minres;r<maxres;r++)
    {
        int first = p->res_fai[r];
        int last = p->res_fai[r+1];
        int at;
        for(at=first;at<last;at++)
        {
            if(!hydrogens && is_hydrogen(p->atoms[at].name))
                continue;


            int hydrophobic = is_hydrophobic(p->atoms[at].name);
            //printf("HPhobic - Atom %s %d is %d\n",p->atoms[at].name,p->atoms[at].resnum,hydrophobic);
            flags[at] = PRESENT;

            if(hphob && hydrophobic)
            {
                //printf("1 Flag for atom %s %d is %d\n",p->atoms[at].name,p->atoms[at].resnum,flags[at]);
                flags[at] = flags[at] | CONTRIB;
                //printf("2 Flag for atom %s %d is %d\n",p->atoms[at].name,p->atoms[at].resnum,flags[at]);
            }

            if(hphil && !hydrophobic)
            {
                //printf("3 Flag for atom %s %d is %d\n",p->atoms[at].name,p->atoms[at].resnum,flags[at]);
                flags[at] = flags[at] | CONTRIB;
                //printf("4 Flag for atom %s %d is %d\n",p->atoms[at].name,p->atoms[at].resnum,flags[at]);
            }

        }
    }//r minres maxres



}// make_asa_list

