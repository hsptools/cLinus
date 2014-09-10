#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"

int bumpcheck(Atom3d *atoms, int *fai, int ires1, int ires2, int jres1, int jres2)
{
    int i, j, ia, ja;

    for (i=ires1; i<ires2; i++)
    {
        int i1 = fai[i];
        int i2 = fai[i+1];
        Atom3d *iCA = &atoms[i1+1];
        
        int jstart;

        if(jres1 > i) jstart = jres1;
        else jstart = i+1;
        
        for (j=jstart; j<jres2; j++)
        {
            int j1 = fai[j];
            int j2 = fai[j+1];

            Atom3d *jCA = &atoms[j1+1];

            if (!close(iCA, jCA, 15.0e0))
                continue;

            for (ia=i1; ia<i2; ia++)
            {
                Atom3d *atom = &atoms[ia];

                for(ja=j1; ja<j2; ja++)
                {
                    Atom3d *jatm = &atoms[ja];
                    //printf("%s %s i %d j %d j-i %d bsep %d\n",atom->name,jatm->name,i,j,j-i,atom->bsep1 + jatm->bsep0);

                    if (close(atom, jatm, -1.0e0))
                    {    
                        //printf("Close\n");

                        if ((j-i) == 1)
                        {
                            if ( !(strcmp(atom->name, " CB ")) && !(strcmp(jatm->name, " H  ")) )
                            {
                                double dist=distance(atom,jatm);
                                double ovrlap = (atom->radius + jatm->radius) - dist;
                                //printf("%lf %lf\n",atom->radius + jatm->radius,dist);
                                if (ovrlap < 0.55)
                                    continue;
                                else
                                    return 1;
                            }
                        } //if ((j-i) == 1)

                        if ((j-i) > 1)
                        {
                            if (!(strcmp(atom->name, " H  ") ||	strcmp(jatm->name, " O  ")) ||
                                !(strcmp(atom->name, " O  ") ||	strcmp(jatm->name, " H  "))) 
                            {
                                double dist=distance(atom,jatm);
                                double ovrlap = (atom->radius + jatm->radius) - dist;
                                //printf("%lf %lf\n",atom->radius + jatm->radius,dist);
                                if (ovrlap < 0.55)
                                   continue;
                                
                            
                            }
                            else
                                return 1;
                           
                        } //if ((j-i) > 1)

                        else if ((atom->bsep1 + jatm->bsep0) > 3)
                        {
                            if (!(strcmp(atom->name, " H  ") || strcmp(jatm->name, " O  ")) ||
                            !(strcmp(atom->name, " H  ") || strcmp(jatm->name, " N  ")) ||
                            !(strcmp(atom->name, " N  ") || strcmp(jatm->name, " H  ")) ||
                            !(strcmp(atom->name, " O  ") || strcmp(jatm->name, " H  "))) 
                            
                            
                            
                            {
                                double dist=distance(atom,jatm);
                                double ovrlap = (atom->radius + jatm->radius) - dist;
                                //printf("%lf %lf\n",atom->radius + jatm->radius,dist);
                                if (ovrlap < 0.55)
                                    continue;
                                
                            }
                            else
                                return 1;
                            
                        }//else if ((atom->bsep1 + jatm->bsep0) > 3)
                    }//close atom jatm
                } // for ja
            }// for ia
        }//for j
    }// for i

    return 0;
}
