#include <stdio.h>
#include <string.h>
#include <math.h>
#include "triangles.h"
#include "linus.h"
#define MAXNEIGHBORS 200
#define M_PI 3.14159265358979323846
#define PRESENT 1
#define CONTRIB 2


  
 /*Return the accessible surface area of a molecule

    Arguments

      o *mol* - instance of linusMol - molecule whose accessible surface
        area is to be calculated

      o *flags* - object - the object returned by the *make_asa_list*
        function in the **construct** module.

      o *probe* - float - the water probe size to use for accessible surface

      o *ndiv* - integer - a number from 1 to 5 representing the number of
        sampling points used in the area calculation.  From least (1) to
        most (5) accurate, this is either 60, 240, 960, 3840, or 15,350
        points.

      o Orig *ext_atoms* - list of coordinates - a list of tuples, each
        tuple containing (x, y, z) of an additional atom to include
        in the ASA calculation
        New - Instead of using list of tuples use a list of coordinates. To acces multiply by 3
        
        

      o *ext_radius* - float - the radius of the atoms in ext_atoms.  
   */
   
   
/* asa_eval calls find_neighbors & tri_buried
 *
 * Given a list of atoms and a Numeric array of flags, for each atom,
 * this function determines the accessible surface area for each atom
 * according to the flags given.  
 
 * If use_data is true, then the ASA for each flagged atom will be stored 
 * in data (a Numeric array if it used).
 
 * If use_ext is nonzero, coords is assumed to hold a list of use_ext
 * atom coordinates to include in the calculation. 
 * ext_rad is the radius of external water atoms in coords.  
 * probe is the probe water radius, 
 * 
 * triangles and ndiv both specify the precision of the calculation.  
 * 
 * To cite this algorithm, use
 * 
 * Shrake, A., and J. A. Rupley.  "Environment and Exposure to Solvent of
 *   Protein Atoms.  Lysozyme and Insulin."  J. Mol. Bio. 79 (1973): 351-
 *   371.
 */

 /* ll_find_neighbors
 *
 * Given an atom index k, this function calculates which atoms are close
 * enough to be relevant in solvent accessible surface area calculations.
 * This set of atom indices are stored in nghlst, and the size of nghlst
 * is returned.  Positive indices in nghlst represent atoms in the protein
 * atom list atoms.  Negative indices (offset by one) represent atom
 * coordinates stored in the external coordinate list.
 */
 
 /* tri_buried
 
 * This function determines whether one of the radial points of an atom
 * is buried by steric contacts with other solvated radii.  tx, ty, tz
 * are points on the solvated radii of the atom in question, i.e. 
 * r = atom radius + solvent radius.  This function checks whether 
 * any of the other solvated atoms collide with this particular point.
 * If it is buried within the molecule, this function returns true.
 */
 
double asa_score(LinusProtein *p, int  *flags, double probe, int ndiv, double *ext_coords, int numext, double ext_radius)
{
    double asa;
    
    if(probe < 0.0) probe = 1.4;
    if(ndiv < 0) ndiv = 3;
    if(ext_radius < 0) ext_radius = 1.4;
    
    Atom3d *atoms = p->atoms;

    //double *ext_coords; 

    int use_data = 0;
    
    double data[p->num_atoms] ;
    
    
    asa = asa_eval(atoms, p->num_atoms, data, ext_coords, flags, probe,
                        ext_radius, use_data, numext, ndiv);
    //printf("Asa is %lf\n",asa);
    return asa;                    

}


double asa_eval(Atom3d * atoms, int numatm, double * data, 
                double *coords, int *flags,
                double probe, double ext_rad,
                int use_data, int numext, int ndiv)
{

     

    //const float triangles[][3];
    double triarea, r, x, y, z, tx, ty, tz, area, atot=0.0;
    int ntrian, i, k, nts, numngh, nghstrt;
    /*
   //printf("Numext %d\n",numext);
    for(i=0;i<numext;i++)
    {
        //printf("%lf ",coords[i]);
    }
    //printf("\n");
    */
    
    int nghlst[MAXNEIGHBORS];
    Atom3d *a;

    switch (ndiv) 
    {
        case 1:
            ntrian = 60;
            break;
        case 2:
            ntrian = 240;
            break;
        case 3:
            ntrian = 960;
            break;
        case 4:
            ntrian = 3840;
            break;
        default:
            ntrian = 15360;
            break;
    }
    //ntrian = 10;
    triarea = 4.0 * M_PI / (double) ntrian;

    for (k=0; k < numatm; k++) 
    {
        //printf("\n\nAtom %s %d\n",atoms[k].name,flags[k]);
        //if (!((flags[k] & CONTRIB) || use_data))  
        if(!(flags[k])) continue;
        
        a = &atoms[k];
       
        numngh = find_neighbors(atoms, coords, flags, probe, ext_rad,k, numatm, numext, nghlst);
    
        /*
         //printf("Atom %i has %i neighbors:\n  ", k, numngh);
            
          for (i=0; i<numngh; i++) 
          {
            printf("%d ", nghlst[i]);
          }
          printf("\n");
        */ 
        
        nghstrt = 0;
        nts = 0;
        r = probe + a->radius;
        //printf("%s %lf ",a->name,a->radius);
        x = a->x; y = a-> y; z = a->z;
    
        
        for (i=0; i<ntrian; i++) 
        {
            tx = ((double) triangles3[i][0])*r + x;
            ty = ((double) triangles3[i][1])*r + y;
            tz = ((double) triangles3[i][2])*r + z;
            //printf("%d %lf %lf %lf %lf %lf %lf\n",i,tx,ty,tz,triangles3[i][0],triangles3[i][1],triangles3[i][2] ); 
            
            if (!tri_buried(atoms, coords,numext, tx, ty, tz, probe, ext_rad,nghlst, numngh, &nghstrt))
            {
                       
                nts++;
            }    

        }
        //printf("%s nts %d\n",atoms[k].name,nts);
        area = (double) nts * triarea * r * r;
        
        if (use_data) 
        {
            data[k] = area;
        }

        if (flags[k] & CONTRIB) 
        {
            atot += area;
        }
        
    }
    //printf("Atot %lf\n",atot);*/
    return atot;
}

int  find_neighbors(Atom3d *atoms, double *coords,int *flags,double probe, 
double ext_rad,int k, int numatm, int numext, int *nghlst)
{
  
    double d, ax, ay, az, bx, by, bz;
    int pos, i, numngh=0;

    Atom3d* a = &atoms[k];
    Atom3d* b;
  
    double arad = a->radius + probe + probe;


    for (i=0; i < k; i++) 
    {
        if (!(flags[i])) continue;

        b = &atoms[i];
        d = distance(a, b);
        //printf("%s %lf %s %lf dist %lf\n",a->name,arad,b->name,b->radius,d);
        if (d < (arad + b->radius)) 
        {
          //printf("Ngh %d NghLst %s %d\n", numngh,b->name,i);
          nghlst[numngh++] = i;
          
        }
    }
    //printf("Before atom %s-  num ngh %d ",a->name,numngh);
    
    for (i=k+1; i < numatm; i++) 
    {
        if (!(flags[i])) continue;
        
        b = &atoms[i];
        d = distance(a, b);
        //printf("%s %lf %s %lf dist %lf\n",a->name,arad,b->name,b->radius,d);
        if (d < (arad + b->radius)) 
        {
            //printf("Ngh %d NghLst %s %d\n", numngh,b->name,i);
            nghlst[numngh++] = i;
        }
    }
    
    //printf("After atom %s- num ngh %d\n",a->name,numngh);
    
    //printf("Numext %d\n",numext); 
    if (!numext) return numngh;

    ax = a->x; ay = a->y; az = a->z;
      
    for (i=0; i < numext/3; i++) 
    {
        int pos = 3*i;
        bx = coords[pos] - ax;
        by = coords[pos+1] - ay;
        bz = coords[pos+2] - az; // Orig bz = coords[pos]   - az

        d = sqrt(bx*bx + by*by + bz*bz);
        if (d < (arad + ext_rad)) 
        {
            nghlst[numngh++] = -(i+1);
        }
    }
    //printf("Extra atom %s-  num ngh %d\n",a->name,numngh);

  return numngh;
}


int tri_buried(Atom3d* atoms, double *coords, int numext, double tx, double ty, double tz, 
double probe, double ext_rad, int *nghlst, int numngh,int *nghstrt)
{
    Atom3d* a;
  
    double dx, dy, dz, d, r;
    int i, k, pos;

    for (k = *nghstrt; k < numngh; k++) 
    {
        //printf("Nghlist %d %d\n",k,nghlst[k]);
        i = nghlst[k];
        
        if (i >= 0) 
        {
            a = &atoms[i];
            dx = a->x - tx; dy = a->y - ty; dz = a->z - tz;
            d = dx*dx + dy*dy + dz*dz;
            r = probe + a->radius;
            //printf("If 1 %s d %lf tx %lf ty %lf tz %lf ax %lf ay %lf az %lf\n",a->name,d,tx,ty,tz,a->x,a->y,a->z);
            if (d < r*r) 
            {
                *nghstrt = k;
                return 1;
            }
        }
        
        else 
        {
           
            pos = 3*(-(i+1));
            dx = coords[pos++] - tx;
            dy = coords[pos++] - ty;
            dz = coords[pos] - tz; // Orig pos
            d = dx*dx + dy*dy + dz*dz;
            r = ext_rad + probe;
            //printf("Else 1 d %lf tx %lf ty %lf tz %lf\n",d,tx,ty,tz);
            if (d < r*r) 
            {
                *nghstrt = k;
                return 1;
            }
        }
        
    }
   
    for (k = 0; k < *nghstrt; k++) 
    {
        i = nghlst[k];
        //printf("Nghlist %d %d\n",k,nghlst[k]);
        
        if (i >= 0) 
        {
            a = &atoms[i];
            dx = a->x - tx; dy = a->y - ty; dz = a->z - tz;
            d = dx*dx + dy*dy + dz*dz;
            r = probe + a->radius;
            //printf("If 2 %s d %lf tx %lf ty %lf tz %lf ax %lf ay %lf az %lf\n",a->name,d,tx,ty,tz,a->x,a->y,a->z);
            if (d < r*r) 
            {
                *nghstrt = k;
                return 1;
            }
        }
      
        else 
        {
            pos = 3*(-(i+1));
            dx = coords[pos++] - tx;
            dy = coords[pos++] - ty;
            dz = coords[pos] - tz;
            d = dx*dx + dy*dy + dz*dz;
            r = ext_rad + probe;
            //printf("Else 1 d %lf tx %lf ty %lf tz %lf\n",d,tx,ty,tz);
            if (d < r*r) 
           {
                *nghstrt = k;
                return 1;
            }
        }

      }
    //printf("No if else\n");  
    return 0;
}


/*
"""Return the accessible surface area of a molecule

    Arguments

      o *filename* - filename or open file - location to which ASA
        information will be written

      o *mol* - instance of linusMol - molecule whose accessible surface
        area is to be calculated

      o *flags* - object - the object returned by the *make_asa_list*
        function in the **construct** module.

      o *probe* - float - the water probe size to use for accessible surface

      o *ndiv* - integer - a number from 1 to 5 representing the number of
        sampling points used in the area calculation.  From least (1) to
        most (5) accurate, this is either 60, 240, 960, 3840, or 15,350
        points.

      o *ext_atoms* - list of coordinates - a list of tuples, each
        tuple containing (x, y, z) of an additional atom to include
        in the ASA calculation

      o *ext_radius* - float - the radius of the atoms in ext_atoms.

      o *title* - string - An optional title for the PDB output
    """

*/



double print_asa_score(LinusProtein *p,int  *flags, double probe, int ndiv, double *ext_coords, int numext, double ext_radius, CntmultList *solv_list,char *filename)
{
    double asa;
    
    int minres,maxres;
    get_res_extents(p,&minres,&maxres);
    
    if(probe < 0.0) probe = 1.4;
    if(ndiv < 0) ndiv = 3;
    if(ext_radius < 0) ext_radius = 1.4;
    
    Atom3d *atoms = p->atoms;

    int use_data = 1;
    
    double data[p->num_atoms] ;
   
    
    
    asa = asa_eval(atoms, p->num_atoms, data, ext_coords, flags, probe,ext_radius, use_data, numext, ndiv);
    printf ("Asa is %lf",asa);
    FILE *file = fopen(filename,"w+");
    
    int i;
    for(i=minres;i<maxres;i++)
    {
        
        double occ = 0.0;
        int start = p->res_fai[i];
        int end = p->res_fai[i+1];
        int j;
        for (j=start;j<end;j++)
        {
            
            if(solv_list[i].is_N ==1)
            {
                if(strcmp(atoms[j].name," N  " ))
                {
                fprintf(file,"ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",j+1,p->atoms[j].name, p->res_names[i],p->res_num[i],p->atoms[j].x,p->atoms[j].y,p->atoms[j].z,(double)solv_list[i].num_N,data[j]);
                }
                else if(strcmp(atoms[j].name," O  " ))
                {
                    fprintf(file,"ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",j+1,p->atoms[j].name, p->res_names[i],p->res_num[i],p->atoms[j].x,p->atoms[j].y,p->atoms[j],(double)solv_list[i].num_O,data[j]);
                }
                else
                {
                    fprintf(file,"ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",j+1,p->atoms[j].name, p->res_names[i],p->res_num[i],p->atoms[j].x,p->atoms[j].y,p->atoms[j].z,occ,data[j]);
                }
        
            }
            else
            {
                fprintf(file,"ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",j+1,p->atoms[j].name, p->res_names[i],p->res_num[i],p->atoms[j].x,p->atoms[j].y,p->atoms[j].z,occ,data[j]);
            }    
        }
      
    }
    
    return asa;    
}
