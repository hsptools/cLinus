#include <stdio.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"


double hbond_score (LinusProtein *p,HbList *hblist, int nhb)
{

    int i, j, nacp;
    double ehb = 0.0e0;
    double dist;

    for(i=0; i<nhb; i++)
    {
        Atom3d donor = p->atoms[hblist[i].donor];
        Atom3d da1 = p->atoms[hblist[i].da1];
        Atom3d da2 = p->atoms[hblist[i].da2];
        nacp = hblist[i].acps_size;

        for (j=0; j<nacp; j++)
        {
            Atom3d acp = p->atoms[hblist[i].acps[j].acp];
            Atom3d acp1 = p->atoms[hblist[i].acps[j].acp1];
            double hbd = hblist[i].acps[j].hbdist;
            double hbm = hblist[i].acps[j].hbdmax;
            double hbs = hblist[i].acps[j].hbene;

            if (!(strcmp(donor.name, " OG ") && strcmp(donor.name, " OG1") &&
            strcmp(donor.name, " OH ") && strcmp(donor.name, " NH1") &&
            strcmp(donor.name, " NZ ") && strcmp(donor.name, " NH2") &&
            strcmp(donor.name, " ND1") && strcmp(donor.name, " ND2") &&
            strcmp(donor.name, " OD2") && strcmp(donor.name, " SG ") &&
            strcmp(donor.name, " NE2") && strcmp(donor.name, " OE2") &&
            strcmp(donor.name, " NE1") && strcmp(donor.name, " OE1") &&
            strcmp(donor.name, " OD1")))
            {
                if (!close(&donor, &acp, hbm))
                  continue;
                if (!(angle(&donor, &acp, &acp1) > 90.0e0))
                  continue;

                dist = distance(&donor, &acp);

                if (dist > hbd)
                  hbs *= (hbm - dist)/(hbm - hbd);

                ehb += hbs;
                /* break; */   /* break was removed to allow for multiple donors
                          and acceptors for one residue. 20040413 */
            }

            else
            {
                if (!close(&donor, &acp, hbm))
                  continue;
                if (!(angle(&donor, &acp, &acp1) > 90.0e0))
                  continue;
                if (!(angle(&da1, &donor, &acp) > 69.0e0))
                  continue;
                if (!(angle(&da2, &donor, &acp) > 69.0e0))
                  continue;
                if (!(fabs(torsion(&acp, &donor, &da1, &da2)) > 130.0e0))
                  continue;

                dist = distance(&donor, &acp);

                if (dist > hbd)
                  hbs *= (hbm - dist)/(hbm - hbd);

                ehb += hbs;
                /* break; */   /* break was removed to allow for multiple donors
                          and acceptors for one residue. 20040413 */
            }
        }
    }

  return -ehb;
}

/*List details of the hydrogen bond score to a file

    Arguments

        o *mol* - instance of LinusMol - molecule whose hydrogen bond
        energy is to be computed

        o *hblist* - object - the object returned by the function
        *make_hbond_list* in the *construct* module

        o *file* - string or file object - name of a file or a file
        previously opened for writing.

    Result

        The hydrogen bond energy as a floating point number
*/


double print_hbond_score(LinusProtein *p, HbList *hblist,int nhb,char *filename)
{
    FILE *file = fopen(filename,"w+");
    fprintf(file,"\n------ Hydrogen Bond Score Listing -----\n");
    int i, j, nacp;
    double ehb = 0.0e0;
    double dist;

    for(i=0; i<nhb; i++)
    {
        Atom3d donor = p->atoms[hblist[i].donor];
        Atom3d da1 = p->atoms[hblist[i].da1];
        Atom3d da2 = p->atoms[hblist[i].da2];
        nacp = hblist[i].acps_size;

        for (j=0; j<nacp; j++)
        {
            Atom3d acp = p->atoms[hblist[i].acps[j].acp];
            Atom3d acp1 = p->atoms[hblist[i].acps[j].acp1];
            double hbd = hblist[i].acps[j].hbdist;
            double hbm = hblist[i].acps[j].hbdmax;
            double hbs = hblist[i].acps[j].hbene;

            if (!(strcmp(donor.name, " OG ") && strcmp(donor.name, " OG1") &&
            strcmp(donor.name, " OH ") && strcmp(donor.name, " NH1") &&
            strcmp(donor.name, " NZ ") && strcmp(donor.name, " NH2") &&
            strcmp(donor.name, " ND1") && strcmp(donor.name, " ND2") &&
            strcmp(donor.name, " OD2") && strcmp(donor.name, " SG ") &&
            strcmp(donor.name, " NE2") && strcmp(donor.name, " OE2") &&
            strcmp(donor.name, " NE1") && strcmp(donor.name, " OE1") &&
            strcmp(donor.name, " OD1")))
            {
                if (!close(&donor, &acp, hbm))
                  continue;
                if (!(angle(&donor, &acp, &acp1) > 90.0e0))
                  continue;

                dist = distance(&donor, &acp);

                if (dist > hbd)
                  hbs *= (hbm - dist)/(hbm - hbd);

                ehb += hbs;

                fprintf(file,"%4d %3s %4s %4d %3s %4s Distance = %lf Score = %lf",
                       donor.resnum,p->res_names[donor.resnum],donor.name,
                       acp.resnum,p->res_names[acp.resnum],acp.name,dist,hbs);
                /* break; */   /* break was removed to allow for multiple donors
                          and acceptors for one residue. 20040413 */
            }

            else
            {
                if (!close(&donor, &acp, hbm))
                  continue;
                if (!(angle(&donor, &acp, &acp1) > 90.0e0))
                  continue;
                if (!(angle(&da1, &donor, &acp) > 69.0e0))
                  continue;
                if (!(angle(&da2, &donor, &acp) > 69.0e0))
                  continue;
                if (!(fabs(torsion(&acp, &donor, &da1, &da2)) > 130.0e0))
                  continue;

                dist = distance(&donor, &acp);

                if (dist > hbd)
                  hbs *= (hbm - dist)/(hbm - hbd);

                ehb += hbs;
                fprintf(file,"%4d %3s %4s %4d %3s %4s Distance = %lf Score = %lf",
                       donor.resnum,p->res_names[donor.resnum],donor.name,
                       acp.resnum,p->res_names[acp.resnum],acp.name,dist,hbs);
                /* break; */   /* break was removed to allow for multiple donors
                          and acceptors for one residue. 20040413 */
            }
        }
    }
    fprintf(file,"\nTotal Hydrogen Bond Score = %lf\n",-ehb);

   fclose(file);

  return -ehb;



}



