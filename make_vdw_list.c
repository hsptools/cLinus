#include <stdio.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"
int make_vdw_list(LinusProtein * p, AtomList * vlist, int size)
{
    int ctr =0;


    int minres,maxres;
    get_res_extents(p,&minres,&maxres);

    double RMIN2SIG = 1.0/(pow(2.0,(1.0/6.0)));
    //printf("RMIN2SIG %lf\n",RMIN2SIG);
    int i;

    for(i=minres; i < maxres-1;i++)
    {
        ////printf("i %d ",i);
        double arad,asigma,brad,bsigma;

        int i1 = p->res_fai[i];
        int i2  = p->res_fai[i+1];

        // do intraresidue first
        int start;
        if(minres)
            start = i1+5;
        else
            start = i1+7;
        //printf("%d %d %d %d\n",minres,maxres,start,i2);    

        //do sidechain with carbonyl oxygen
        // biggest sidechain 16
        if(start<i2)
        {
            Atom3d atom = p->atoms[i1+3];
            //printf("%s %d\n",atom.name,atom.resnum);
            get_atom_vdw_radius(p->res_names[i],atom.name,&arad);
            asigma = (arad * RMIN2SIG);
            int ja;
            for(ja = start;ja<i2;ja++)
            {
                Atom3d jatom = p->atoms[ja];
                //printf("%d %d %d %s %d\n",start,i2,ja,jatom.name,jatom.resnum);
                if(!strcmp(jatom.name, " H  ") || !strcmp(jatom.name, " HA "))
                    continue;
                get_atom_vdw_radius(p->res_names[i],jatom.name,&brad);
                bsigma = (brad*RMIN2SIG);

                AtomList vlist1;//={atom,jatom,asigma+bsigma};
                vlist1.a1 = i1+3;
                vlist1.a2 = ja;
                
                vlist1.value = asigma+bsigma;
                if(size)
                {
                    vlist[ctr] = vlist1;
                }
                ctr++;
            }
        }

        int j;
        for (j=i+1;j<maxres;j++)
        {
            int j1 = p->res_fai[j];
            int j2 = p->res_fai[j+1];
            int ia;

            for(ia=i1;ia<i2;ia++)
            {
                Atom3d atom1 = p->atoms[ia];
                if(!strcmp(atom1.name, " H  ") || !strcmp(atom1.name, " HA "))
                    continue;
                get_atom_vdw_radius(p->res_names[i],atom1.name,&arad);
                asigma = (arad * RMIN2SIG);


                int ja;
                for(ja=j1;ja<j2;ja++)
                {
                    Atom3d jatom1 = p->atoms[ja];
                    if(!strcmp(jatom1.name, " H  ") || !strcmp(jatom1.name, " HA "))
                        continue;
                    get_atom_vdw_radius(p->res_names[j],jatom1.name,&brad);
                    bsigma = (brad*RMIN2SIG);

                    if(j-i > 1)
                    {
                        ////printf("if j-i>1\n");
                        AtomList vlist1;
                                             
                        vlist1.a1 = ia;
                        vlist1.a2 = ja;
                        vlist1.value = asigma+bsigma;
                        if(size)
                        {
                            vlist[ctr] = vlist1;
                        }
                        ctr++;;
                       

                    }

                    else if(atom1.bsep1 + jatom1.bsep0 > 3)
                    {

                        AtomList vlist1;
                        vlist1.a1 = ia;
                        vlist1.a2 = ja;
                        vlist1.value = asigma+bsigma;
                        
                        if(size)
                        {
                            vlist[ctr] = vlist1;
                        }

                        ctr++;
                    }

                    else if(!strcmp(atom1.name, " C  ") && !strcmp(jatom1.name, " O  "))
                    {

                        ////printf("atom name ja\n");
                        AtomList vlist1;
                        vlist1.a1 = ia;
                        vlist1.a2 = ja;
                        vlist1.value = asigma+bsigma;
                        if(size)
                        {
                            vlist[ctr] = vlist1;
                        }

                        ctr++;;

                    }
                }

            }
        }
    }
    return ctr;
}

