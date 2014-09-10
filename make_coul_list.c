#include <stdio.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"

// makes list of atom1,atom2,q*q
// for NH and CO in same residue
// q values are last record of linusRes.py

int make_coul_list(LinusProtein * p, AtomList * qlist, int size)
{
    int ctr =0;

    int minres,maxres;
    get_res_extents(p,&minres,&maxres);
    //printf("Min Max %d %d\n",minres,maxres);
    double achrg,bchrg;

    int i;
    for(i=minres;i < (maxres-1);i++)
    {
        int i1 = p->res_fai[i];
        int i2 = p->res_fai[i+1];
        //printf("i %d\n",i);
        //int j1 = p->res_fai[i+1];
        //int j2 = p->res_fai[i+2];

        //int k1 = p->res_fai[i-1];
        //int k2 = p->res_fai[i];

        int ia;
        for(ia = i1; ia <i2;ia++)
        {
            Atom3d atom = p->atoms[ia];

            //printf("%d %d Atom1 name %s %d \n",i,ia,atom.name,atom.resnum,p->res_names[i]);
            if( !strcmp(atom.name," N  ") )
            {

                get_atom_chrg(p->res_names[i],atom.name,&achrg);

                int ja;
                for(ja=i1;ja<i2;ja++)
                {
                    Atom3d jatom = p->atoms[ja];

                    //printf("%d Atom2 name %s %d\n",ja,jatom.name,jatom.resnum,p->res_names[i]);
                    if( !strcmp(jatom.name," C  ") )
                    {
                        //printf("NC\n");
                        get_atom_chrg(p->res_names[i],jatom.name,&bchrg);
                        AtomList qlist1;
                        qlist1.a1 = ia;
                        qlist1.a2 = ja;
                        qlist1.value = achrg*bchrg;
                        if(size) qlist[ctr] = qlist1;
                        ctr++;
                    }
                }
            }

            if( !strcmp(atom.name," N  ") )
            {
                get_atom_chrg(p->res_names[i],atom.name,&achrg);

                int ja;
                for(ja=i1;ja<i2;ja++)
                {
                    Atom3d jatom = p->atoms[ja];
                    if( !strcmp(jatom.name," O  ") )
                    {
                        //printf("NO\n");
                        get_atom_chrg(p->res_names[i],jatom.name,&bchrg);
                        AtomList qlist1;
                        qlist1.a1 = ia;
                        qlist1.a2 = ja;
                        qlist1.value = achrg*bchrg;
                        if(size) qlist[ctr] = qlist1;
                        ctr++;
                    }
                }
            }


            if( !strcmp(atom.name," H  ") )
            {
                get_atom_chrg(p->res_names[i],atom.name,&achrg);

                int ja;
                for(ja=i1;ja<i2;ja++)
                {
                    Atom3d jatom = p->atoms[ja];
                    if( !strcmp(jatom.name," C  ") )
                    {
                        //printf("HC\n");
                        get_atom_chrg(p->res_names[i],jatom.name,&bchrg);
                        AtomList qlist1;
                        qlist1.a1 = ia;
                        qlist1.a2 = ja;
                        qlist1.value = achrg*bchrg;
                        if(size) qlist[ctr] = qlist1;
                        ctr++;
                    }
                }
            }

            if( !strcmp(atom.name," H  ") )
            {
                get_atom_chrg(p->res_names[i],atom.name,&achrg);

                int ja;
                for(ja=i1;ja<i2;ja++)
                {
                    Atom3d jatom = p->atoms[ja];
                    if( !strcmp(jatom.name," O  ") )
                    {
                        //printf("HO\n");
                        get_atom_chrg(p->res_names[i],jatom.name,&bchrg);
                        AtomList qlist1;
                        qlist1.a1 = ia;
                        qlist1.a2 = ja;
                        qlist1.value = achrg*bchrg;
                        if(size) qlist[ctr] = qlist1;
                        ctr++;
                    }
                }
            }



        }

    }



return ctr;
}
