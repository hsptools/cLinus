#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"
#include "LinusRes.h"
#include "Chidict.h"

// Does not read generic pdb file. Reads only a specific file in tetra ala format


void protein_from_pdb(char * file_name, LinusProtein *p)
{
    Atom3d * atoms;                                  // array of Atoms

    char atm[5];
    char aname[6];
    char rname[4];                                  //replaced res_name
    int anum;                                       //replaced serial
    int rnum;                                       //replaced res_seq
    int num_atoms=0;                                // number of lines in .pdb file //replaced count
    int num_res = 0;
    double x;
    double y;
    double z;
    double occupancy;
    double tempFactor;
    char buff[BUFSIZ];

    char *aname_seen;                               // atom names seen array
    //int res_start_end[2];                          // array containing start/end info for residue in .pdb //replaced numres
    int start, end, res_size;
    int i, j; // loop variables
    int current_res = -1;
    int na = 0; // atom indices
    int nr = -1; // residue indices
    int *fai; // na array // replaced f
    int *rnum_arr; // replaced respdbnum

    strcpy(p->name, "PROT"); // name protein PROT
    strcpy(p->filename, file_name); // protein filename is first command line argument

    //printf("Name %s, Filename %s\n", p->name, p->filename);

    FILE *file = fopen(p->filename, "r");

    // Get file line count to allocate matrices
    p->num_chi_atoms=0;
    while (fgets(buff, sizeof buff, file) != NULL)
    {
        if (sscanf(buff, "%s %d %s %s %d %lf %lf %lf %lf %lf \n",
         atm, &anum, aname, rname, &rnum, &x, &y, &z,
         &occupancy, &tempFactor) == 10)
        {
            //printf("%s\n",rname);
            if(rnum != current_res)
            {
                num_res++;
                num_atoms+=get_res_size(rname);
                //printf("Res %s res size %d\n", rname, get_res_size(rname));
                p->num_chi_atoms += get_chi_res_size(rname);
                current_res = rnum;
            }
            
        }
    }

    //printf("Num of chi atoms %d\n", p->num_chi_atoms);
    //printf("Num of atoms is %d & Num_residues is %d\n", num_atoms,num_res);

    // Dynamic allocations based on # lines in .pdb file & residue information

    p->res_names = (char **)malloc(num_res * sizeof(char *));
    p->mesostates = (char **)malloc(num_res * sizeof(char *));
    for (i = 0; i < num_res; i++)
    {
        p->res_names[i] = (char *)malloc(sizeof(char)*4);
        p->mesostates[i] = (char *)malloc(sizeof(char)*3);
        strcpy(p->mesostates[i],"--");
    }

    p->atoms = malloc(num_atoms * sizeof(Atom3d));
  
   
   
   
    //printf("Start phi psi ome chi allocation\n");
    p->phi_atoms = malloc(num_res * sizeof(int));
    p->psi_atoms = malloc(num_res * sizeof(int));
    p->omega_atoms = malloc(num_res * sizeof(int));
       
    //printf("Finished phi psi \nStart chi allocation\n");
    
    p->chi_atoms = malloc(p->num_chi_atoms * sizeof(int));
        
    //printf("Success allocating phi psi ome chi\n");
    
    p->res_num = malloc(num_atoms * sizeof(int));

    p->res_fai = malloc((num_res+1) * sizeof(int));

    p->res_helix_wt = malloc((num_res) * sizeof(double));
    p->res_strand_wt = malloc((num_res) * sizeof(double));
    p->res_turn1_wt = malloc((num_res) * sizeof(double));
    p->res_turn2_wt = malloc((num_res) * sizeof(double));
    p->res_coil_wt = malloc((num_res) * sizeof(double));
    p->res_PII_wt = malloc((num_res) * sizeof(double));
    p->orig_fp_distance = malloc((num_res) * sizeof(double));
    p->orig_sp_angle = malloc((num_res) * sizeof(double));
    p->orig_tp_torsion = malloc((num_res) * sizeof(double));
    p->res_f_mean = malloc((num_res) * sizeof(double));
    p->res_f_sd = malloc((num_res) * sizeof(double));
    p->res_y_mean = malloc((num_res) * sizeof(double));
    p->res_y_sd = malloc((num_res) * sizeof(double));
    p->res_acprat = malloc((num_res) * sizeof(double));

    p->pmeso_wt = (double **) malloc((num_res-1) * sizeof(double *));
    for(i=0;i<num_res-1;i++)
    {
        p->pmeso_wt[i] = (double *) malloc(144 * sizeof(double));
    }

    //printf("Finished allocating memory\n");


    // Fill up LinusProtein
    i = 0;
    rewind(file);
    while (fgets(buff, sizeof buff, file) != NULL)
    {
        if (sscanf(buff, "%s %d%*c%4[ 0-9A-Z]%*c%s %d %lf %lf %lf %lf %lf \n",
            atm, &anum, aname, rname, &rnum, &x, &y, &z,
            &occupancy, &tempFactor) == 10)
        {
            if (!strcmp(atm, "ATOM"))
            {
                //printf("-%s-%d-%s-%s-%d-%lf-%lf-%lf-%lf-%lf-\n", atm,anum,aname,rname,rnum,x,y, z,occupancy,tempFactor);
                //na = na + 1;
                if (rnum != current_res) // handles atoms up to max(rnum)-1
                {
                   nr = nr + 1;
                   j = nr;
                   p->res_fai[j] = na;
                   na += get_res_size(rname);
                   strcpy(p->res_names[j], rname);
                   p->res_num[j] = rnum;
                   current_res = rnum;

                   // New use of start/end
                   //start = na;
                   //end = get_res_size(rname); // this may be case specific
                }

                if (!strcmp(aname, "OXT"))
                {
                   na = na - 1;
                   continue;
                }

              
               
                i++;

            } // end atm if-loop

        } // end sscanf if-loop

    } // end while loop
    //fclose(file);
    //printf("Checking 2\n");
    p->num_atoms = num_atoms;;
    p->num_res = num_res;
    p->res_fai[num_res] = num_atoms;
    
    
    int atm_ctr=0;
    int s=0;
    for (i=0;i<num_res;i++)
    {
        int start, end;
        get_res_range(p->res_names[i],&start,&end);
        for(j=start;j<end;j++)
        {
            strcpy(p->atoms[atm_ctr].name,atm_name_arr[j]);
            p->atoms[atm_ctr].resnum = i;
            
            //printf("%d %s %s %d\n",atm_ctr,p->atoms[atm_ctr].name,p->res_names[i],i);
            atm_ctr++;
        }
    
    }
    
    j=-1;
    rewind(file);
    while (fgets(buff, sizeof buff, file) != NULL)
    {
        if (sscanf(buff, "%s %d%*c%4[ 0-9A-Z]%*c%s %d %lf %lf %lf %lf %lf \n",
            atm, &anum, aname, rname, &rnum, &x, &y, &z,
            &occupancy, &tempFactor) == 10)
        {
            if (!strcmp(atm, "ATOM"))
            {
                //printf("-%s-%d-%s-%s-%d-%lf-%lf-%lf-%lf-%lf-\n", atm,anum,aname,rname,rnum,x,y, z,occupancy,tempFactor);
                if (rnum != current_res) // handles atoms up to max(rnum)-1
                {
                   j = j + 1;
                   current_res = rnum;
                   
                }   
                int idx = get_linusRes_index(rname,aname); // get index for linusRes arrays

                if(idx == -1)
                {
                    //printf("Error unknown atom name!");
                    exit(1);
                }

                //strcpy(p->atoms[i].name,aname);
                
                int atm_ind = get_atom_ind_with_name(p->atoms,p->res_fai[j],p->res_fai[j+1],aname);
                //printf("%d %s pdb %s atoms %s fai %d fai+1 %d atm_ind %d\n",j,p->res_names[j],aname,p->atoms[atm_ind].name,p->res_fai[j],p->res_fai[j+1],atm_ind);
                p->atoms[atm_ind].radius = vdw_radius_arr[idx];
                p->atoms[atm_ind].bsep0 = bsep0_arr[idx];
                p->atoms[atm_ind].bsep1 = bsep1_arr[idx];
                strcpy(p->atoms[atm_ind].first_name,first_name_arr[idx]);
                p->atoms[atm_ind].fp = first_offset_arr[idx];
                strcpy(p->atoms[atm_ind].second_name,second_name_arr[idx]);
                p->atoms[atm_ind].sp = second_offset_arr[idx];
                strcpy(p->atoms[atm_ind].third_name, third_name_arr[idx]);
                p->atoms[atm_ind].tp = third_offset_arr[idx];
                p->atoms[atm_ind].sasa = 0.0;
                p->atoms[atm_ind].sasarad = 0.0;
                //printf("First %s Second %s Third %s\n", first_name_arr[idx], second_name_arr[idx],third_name_arr[idx]);
                //printf("Res %s Resnum %d Atom %s Atomnum %d First %s Second %s Third %s\n",rname,rnum, aname, i,p->atoms[i].first_name, p->atoms[i].second_name, p->atoms[i].third_name);

                p->atoms[atm_ind].x = x;
                p->atoms[atm_ind].y = y;
                p->atoms[atm_ind].z = z;
                

                //printf("X %lf Y %lf Z %lf\n",p->atoms[atm_ind].x,p->atoms[atm_ind].y,p->atoms[atm_ind].z);
            } // end atm if-loop

        } // end sscanf if-loop

    } // end while loop
    fclose(file);
     // Finish filling atoms array info
    
   
   
   
   
   
   
   
   
    //From LinusMol init

    int r;
    for(r=0;r<num_res;r++)
    {
        p->res_helix_wt[r] = 0.1;
        p->res_strand_wt[r] = 0.4;
        p->res_turn1_wt[r] = 0.05;
        p->res_turn2_wt[r] = 0.05;
        p->res_coil_wt[r] = 0.20;
        p->res_PII_wt[r] = 0.10;
        p->res_f_mean[r] = 0.0;
        p->res_f_sd[r] = 99.0;
        p->res_y_mean[r] = 0.0;
        p->res_y_sd[r] = 99.0;
        p->res_acprat[r] = 0.0;

    }

    for(r=0;r<144;r++)
    {
        p->pmeso_grid[r] = 0.0;
    }

    p->pmeso_grid[35] = 1.0; //set Ck as probability 1.0

    int q;
    for(r=0;r<num_res-2;r++) // subtract 2 because ACE and NME don't have pmeso_weights
    {
        for(q=0;q<144;q++)
        {
            p->pmeso_wt[r][q] = p->pmeso_grid[q];
        }
    }


   /*
    int g;
    for (g=0;g<num_res;g++)
    {
        printf("Testing residues %s\n",p->res_names[g]);
    }
    */
    

    // Create and initilaize Z matrix
    int n;
    for (n=0;n < num_atoms;n++)
    {
       // printf("Atom name %s Fp %d  %s Resnum %d resnum %d\n",p->atoms[n].name, p->atoms[n].fp,p->atoms[n].first_name, p->atoms[n].resnum, num_res);
        //printf("Atom name %s %d\n",p->atoms[n].name, p->atoms[n].resnum);
        int fp = p->atoms[n].fp + p->atoms[n].resnum;
        if(fp >= 0 && fp < num_res)
        {
            int m = get_atom_ind_with_name(p->atoms, p->res_fai[fp],p->res_fai[fp+1],p->atoms[n].first_name);
            p->atoms[n].fp_distance = distance(&p->atoms[n],&p->atoms[m]);
            p->atoms[n].fp_ind = m;
            //printf(" 1 parent %d %s %lf\n",p->atoms[n].fp_ind,p->atoms[p->atoms[n].fp_ind].name,p->atoms[n].fp_distance);
        }
        else
        {
            p->atoms[n].fp_ind = -1;
            p->atoms[n].fp_distance = 0.0;
        }
        
       
        //printf("1 second_parent %s\n",p->atoms[n].second_parent->name);
        //printf("1 third_parent %s\n",p->atoms[n].third_parent->name);

        int sp = p->atoms[n].sp + p->atoms[n].resnum;
        if(sp >= 0 && sp < num_res)
        {
             //atoms[n].sp_ind = get_atom_ind_with_name(p->atoms, p->res_fai[sp],p->res_fai[sp+1],p->atoms[n].second_name);
             int m = get_atom_ind_with_name(p->atoms, p->res_fai[sp],p->res_fai[sp+1],p->atoms[n].second_name);
             p->atoms[n].sp_angle = angle(&p->atoms[n],&p->atoms[p->atoms[n].fp_ind],&p->atoms[m]);
             p->atoms[n].sp_ind = m;
             //printf(" 2 parent %s %lf\n",p->atoms[m].name,p->atoms[n].sp_angle);
        }

        else
        {
            p->atoms[n].sp_ind = -1;
            p->atoms[n].sp_angle = 0.0;
        }
       
        
        
        int tp = p->atoms[n].tp + p->atoms[n].resnum;
        if(tp >= 0 && tp < num_res)
        {
            
            p->atoms[n].tp_ind = get_atom_ind_with_name(p->atoms, p->res_fai[tp],p->res_fai[tp+1],p->atoms[n].third_name);
            p->atoms[n].tp_torsion = torsion(&p->atoms[n],&p->atoms[p->atoms[n].fp_ind],&p->atoms[p->atoms[n].sp_ind],&p->atoms[p->atoms[n].tp_ind]);
            //printf(" 3 parent %s %lf\n",p->atoms[p->atoms[n].tp_ind].name,p->atoms[n].tp_torsion);
        }
        else
        {
            p->atoms[n].tp_ind = -1;
            p->atoms[n].tp_torsion = 0.0;
        }
      
           

    }
    /*
    
       for(i=0;i<p->num_atoms;i++)
    {
            printf(" Atom %s fp %s sp %s tp %s\n",p->atoms[i].name, p->atoms[p->atoms[i].fp_ind].name,
            p->atoms[p->atoms[i].sp_ind].name,p->atoms[p->atoms[i].tp_ind].name);
          

    }
    
*/
    commit_coords_range(p->atoms, 0,num_atoms);

    // Fill up Phiatoms, Psiatoms, Omeatoms & Chiatoms


    int k;

    int chi_num=0;


    for(k=0;k<num_res;k++)
    {
        p->phi_atoms[k] = get_atom_ind_with_name(p->atoms, p->res_fai[k],p->res_fai[k+1]," C  ");
        p->psi_atoms[k] = get_atom_ind_with_name(p->atoms, p->res_fai[k],p->res_fai[k+1]," O  ");
        p->omega_atoms[k] = get_atom_ind_with_name(p->atoms, p->res_fai[k],p->res_fai[k+1]," CA ");
        
        
        //printf("Phi atom name %s res %d %lf %lf %lf tor %lf\n",p->atoms[p->phi_atoms[k]].name, p->atoms[p->phi_atoms[k]].resnum,p->atoms[p->phi_atoms[k]].x,p->atoms[p->phi_atoms[k]].y,p->atoms[p->phi_atoms[k]].z,p->atoms[p->phi_atoms[k]].tp_torsion);
       // printf("psi atom name %s res %d %lf %lf %lf tor %lf\n",p->atoms[p->psi_atoms[k]].name, p->atoms[p->psi_atoms[k]].resnum,p->atoms[p->psi_atoms[k]].x,p->atoms[p->psi_atoms[k]].y,p->atoms[p->psi_atoms[k]].z,p->atoms[p->psi_atoms[k]].tp_torsion);
        //printf("omega atom name %s res %d %lf %lf %lf tor %lf\n",p->atoms[p->omega_atoms[k]].name, p->atoms[p->omega_atoms[k]].resnum,p->atoms[p->omega_atoms[k]].x,p->atoms[p->omega_atoms[k]].y,p->atoms[p->omega_atoms[k]].z,p->atoms[p->omega_atoms[k]].tp_torsion);
      

        // Fill up chi atoms

        int chi_idx = get_chi_res_index(p->res_names[k]);
        //printf("Chi idx %d\n",chi_idx);
       // printf("Residue %s %d\n",p->res_names[k],k);
        //printf("Residue %s %d\n",p->res_names[k],k);
        if(chi_idx != -1)
        {
            int u;
            int num_chi = chi_size[chi_idx]/2;
            //printf("Num chi %d\n", num_chi);
            for(u=0;u<num_chi;u++)
            {
                //get_atom_with_name(p->chi_atoms[chi_num],p->atoms,p->res_fai[k],p->res_fai[k+1],chi_atoms[chi_idx][u]);
                int m = get_atom_ind_with_name(p->atoms,p->res_fai[k],p->res_fai[k+1],chi_atoms[chi_idx][u]);
                
                //printf("Chi atoms for res %s : %s\n", p->res_names[k],p->atoms[p->chi_atoms[chi_num]].name);
                if(m > -1)
                {
                    p->chi_atoms[chi_num] = m;
                    chi_num++;
                }

            }

        }
        

    }
    p->num_chi_atoms = chi_num;

 
    /*

   free(atm);
   free(aname);
   free(rname);
   free(p->atoms);
   free(p->res_names);
   free(p->res_fai);
   for (i = 0; i < num_res; i++)
   {
       free(p->res_names[i]);
   }
   free(p->res_num);



*/

}
