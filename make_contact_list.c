#include <stdio.h>
//#include "Atom3d.h"
#include "linus.h"
int make_contact_list(LinusProtein * p, int wmin, int wmax, double cscale, AtomList * clist, int size)
{
    if(cscale < 0.0)
        cscale = 1.0;
    int ctr =0;

    int minres,maxres;
    get_res_extents(p,&minres,&maxres);

    double econ = 0.1*cscale;

    int i;
    for(i=minres; i < maxres-1;i++)
    {
        //printf("res %s ",p->res_names[i]);
        int start1 = p->res_fai[i];
        int end1  = p->res_fai[i+1];

        int m;
        for(m=start1;m<end1;m++)
        {
            Atom3d atom = p->atoms[m];
            //printf("Atom %s ",atom.name);
            double arad;
            get_atom_contact_radius(p->res_names[i],atom.name,&arad);
            //printf(" arad %lf \n",arad);
            if(arad <=0.000001)
                continue;

            int start2 = i + wmin;
            int end2;
            if(maxres-1 > i+wmax)
                end2 = i+wmax;
            else
                end2 = maxres-1;
            //printf("start2 end 2 %d %d",start2,end2);

            int j;
            for(j=start2;j<end2;j++)
            {
                int jstart = p->res_fai[j];
                int jend = p->res_fai[j+1];
                int n;
                for(n=jstart;n<jend;n++)
                {
                    Atom3d jatom = p->atoms[n];
                    //printf("jatom %s ",jatom.name);
                    double brad;
                    get_atom_contact_radius(p->res_names[j],jatom.name,&brad);
                    //printf(" brad %lf \n",brad);
                    if(brad >=0.000001)
                    {
                        //printf("Atom %s %d Jatom %s %d Value %lf\n",atom.name,atom.resnum,jatom.name,jatom.resnum,arad+brad);
                        AtomList clist1;
                        clist1.a1 = m;
                        clist1.a2 = n;
                        clist1.value = arad+brad;
                        clist1.econ = econ;
                        if(size) clist[ctr] = clist1;
                        ctr++;
                    }

                }

            }

        }
    }
    return ctr;

}
