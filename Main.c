#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include "linus.h"
#include "linusSecStr.h"

int main(int argc, char** argv) 
{
  create_poly_ala_tet_hexo(argv[1],argv[2]);
 //decorate_box(argv[1],argv[2],argv[3]);
//create_hexo_mesoshort(argv[1],argv[2]);


  /*
    double s[6];
    gen_seeds(s);

    LinusProtein *p = malloc(sizeof(LinusProtein));
    protein_from_pdb(argv[1],p);
    
    p->atoms[p->phi_atoms[2]].tp_torsion = 172.347;
    p->atoms[p->phi_atoms[2]].sp_angle = uniform(111.0,113.0,s);
    p->atoms[p->psi_atoms[2]].tp_torsion = set_sy(-155.611);
    p->atoms[p->omega_atoms[2]].tp_torsion = 182.955;
    
    p->atoms[p->phi_atoms[3]].tp_torsion = 179.887;
    p->atoms[p->phi_atoms[3]].sp_angle = uniform(111.0,113.0,s);
    p->atoms[p->psi_atoms[3]].tp_torsion = set_sy(-161.573);
    p->atoms[p->omega_atoms[3]].tp_torsion = 182.766;
    
    p->atoms[p->phi_atoms[4]].tp_torsion = 172.972;
    p->atoms[p->phi_atoms[4]].sp_angle = uniform(111.0,113.0,s);
    p->atoms[p->psi_atoms[4]].tp_torsion = set_sy(-165.714);
    p->atoms[p->omega_atoms[4]].tp_torsion = 183.115;
    
    p->atoms[p->phi_atoms[5]].tp_torsion = -167.524;
    p->atoms[p->phi_atoms[5]].sp_angle = uniform(111.0,113.0,s);
    p->atoms[p->psi_atoms[5]].tp_torsion = set_sy(-165.420);
    p->atoms[p->omega_atoms[5]].tp_torsion = 175.542;
     
    ztox(p->atoms, 0,p->num_atoms);
    
    int wmin = 2; //2
    int wmax = p->num_res - 2;
    Hbdpar hbparms;
    default_hbdpar(&hbparms);
    hbparms.use_hbond = 1;
    hbparms.use_sidechain_hbond = 0;
   
   
   // Make Hbond list
   
    HbList *hblist;
    int max_acp_size=0,size=0;
    int num_hb = make_hbond_list(p,hblist,wmin,wmax,hbparms.hbond_distance,hbparms.hbond_probe,hbparms.hbond_torsion,
               hbparms.hbond_score_short,hbparms.hbond_score_long,hbparms.sidechain_hbond_distance,hbparms.sidechain_hbond_torsion,
               hbparms.sidechain_hbond_score,hbparms.use_hbond,hbparms.use_sidechain_hbond, &size,&max_acp_size);

    hblist = malloc(num_hb*sizeof(HbList));
    int i;   
    for(i=0;i<size;i++)
    {
        hblist[i].acps = malloc(max_acp_size*sizeof(Acps));
    }   
    make_hbond_list(p,hblist,wmin,wmax,hbparms.hbond_distance,hbparms.hbond_probe,hbparms.hbond_torsion,
       hbparms.hbond_score_short,hbparms.hbond_score_long,hbparms.sidechain_hbond_distance,hbparms.sidechain_hbond_torsion,
       hbparms.sidechain_hbond_score,hbparms.use_hbond,hbparms.use_sidechain_hbond,&num_hb,&max_acp_size);
    
    
    if (!bumpcheck(p->atoms, p->res_fai, 1, p->num_res-2, 2, p->num_res-1)) 
    {
        printf("phitp1 %lf phisp1 %lf psi1 %lf ome1 %lf ",p->atoms[p->phi_atoms[2]].tp_torsion,p->atoms[p->phi_atoms[2]].sp_angle,p->atoms[p->psi_atoms[2]].tp_torsion,p->atoms[p->omega_atoms[2]].tp_torsion);
        printf("phitp2 %lf phisp2 %lf psi2 %lf ome2 %lf ",p->atoms[p->phi_atoms[3]].tp_torsion,p->atoms[p->phi_atoms[3]].sp_angle,p->atoms[p->psi_atoms[3]].tp_torsion,p->atoms[p->omega_atoms[3]].tp_torsion);
        printf("phitp3 %lf phisp3 %lf psi3 %lf ome3 %lf ",p->atoms[p->phi_atoms[4]].tp_torsion,p->atoms[p->phi_atoms[4]].sp_angle,p->atoms[p->psi_atoms[4]].tp_torsion,p->atoms[p->omega_atoms[4]].tp_torsion);
        printf("phitp4 %lf phisp4 %lf psi4 %lf ome4 %lf ",p->atoms[p->phi_atoms[5]].tp_torsion,p->atoms[p->phi_atoms[5]].sp_angle,p->atoms[p->psi_atoms[5]].tp_torsion,p->atoms[p->omega_atoms[5]].tp_torsion);
                           
        FILE *file;
        int numint,numvirt,bbtot,numbb;
        double *wat_list;
        int ctr_wat_list =0;
        CntmultList cntmult_list[p->num_res];
        int wlist = mk_virt_bb_cntmult_loosHbd(p,hblist,num_hb,0,0,file,0,1.25,&numint,&numvirt,&bbtot,&numbb,wat_list,0,cntmult_list);
        wat_list = (double *) malloc(wlist*sizeof(double));
        mk_virt_bb_cntmult_loosHbd(p,hblist,num_hb,0,0,file,0,1.25,&numint,&numvirt,&bbtot,&numbb,wat_list,wlist,cntmult_list);
        //printf("Finished mk virt\n");
      
        int hbs = numbb;
        //printf("hbs %d ",hbs);
        if(hbs > 7)
        {
            printf("Hbonds %d\n",hbs);
            
        }
    }           

    

    for (i=0;i<p->num_atoms;i++)
    {
        printf("%d %s %lf %lf %lf %lf %lf %lf\n",i,p->atoms[i].name,p->atoms[i].x,p->atoms[i].y,p->atoms[i].z,p->atoms[i].fp_distance,p->atoms[i].sp_angle,p->atoms[i].tp_torsion);
    }
  
        int i;
    printf("Waters\n");
    int num_waters = add_waters(p);
    
    /*
    for(i=0;i<num_waters;i++)
    {
        printf("Res %s Fp %s Sp %s Tp %s Dist %lf Ang %lf Tor %lf X %lf Y %lf Z %lf\n",p->res_names[p->waters[i].resnum],p->atoms[p->waters[i].fp_ind].name, p->atoms[p->waters[i].sp_ind].name,p->atoms[p->waters[i].tp_ind].name,p->waters[i].fp_distance, p->waters[i].sp_angle, p->waters[i].tp_torsion, p->waters[i].x, p->waters[i].y, p->waters[i].z);
    }
*/
    /*find_native_like(atof(argv[1]),atof(argv[2]),atof(argv[3]),atof(argv[4]),
                     atof(argv[5]),atof(argv[6]),atof(argv[7]),atof(argv[8]),
                     atof(argv[9]),atof(argv[10]),atof(argv[11]),atof(argv[12]),
                     atof(argv[13]),argv[14],argv[15]); 


// ./findnative -71.8 -40.9 179.9 -54.4 -49.1 173.1 -60.8 -41.3 179.7 -66.7 -33.3 179.7 30 test10k_successful_decorations vvme

                     
    //decorate_box(argv[1],argv[2],argv[3]);
    //decorate_box(argv[1],argv[2],argv[3]);

    //create_poly_ala_tet(argv[1],argv[2]);

/*
    LinusProtein *prot = malloc(sizeof(LinusProtein));
   
    protein_from_pdb(argv[1],prot);

    double DH1,DH2,DH3;
    int a[6];
    int m;
    for(m=0;m < prot->num_res;m++)
    {
        a[m] = get_atom_ind_with_name(prot->atoms,prot->res_fai[m],prot->res_fai[m+1]," CA ");
    }
    DH1 = _tor(&prot->atoms[a[0]],&prot->atoms[a[1]],&prot->atoms[a[2]],&prot->atoms[a[3]]);
    DH2 = _tor(&prot->atoms[a[1]],&prot->atoms[a[2]],&prot->atoms[a[3]],&prot->atoms[a[4]]);
    DH3 = _tor(&prot->atoms[a[2]],&prot->atoms[a[3]],&prot->atoms[a[4]],&prot->atoms[a[5]]);

    printf("%lf %lf %lf\n",DH1,DH2,DH3);


    //

   /* LinusProtein *prot = malloc(sizeof(LinusProtein));
   
    protein_from_pdb(argv[1],prot);

    
     HbList *hblist;
    int i,j,size=0,max_acp_size=0;
       
    int num_hb = make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&size,&max_acp_size);
    hblist=malloc(size*sizeof(HbList));
   
    for(i=0;i<size;i++)
    {
        hblist[i].acps = malloc(max_acp_size*sizeof(Acps));
    }

    make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&num_hb,&max_acp_size);
    
    
    int numint_loos, numvirt_loos, bbtot,numbb;
    CntmultList solv_list[prot->num_res];
    double *wat_list;
    
    FILE *file;
    
    int wlist = mk_virt_bb_cntmult_loosHbd(prot,hblist,size,0,0,file,0,1.25,&numint_loos,&numvirt_loos,&bbtot,&numbb,wat_list,0,solv_list);
    wat_list = (double *) malloc(wlist*sizeof(double));
    numint_loos= numvirt_loos=bbtot=numbb=0;
    mk_virt_bb_cntmult_loosHbd(prot,hblist,size,0,0,file,0,1.25,&numint_loos,&numvirt_loos,&bbtot,&numbb,wat_list,wlist,solv_list);
    

    
    int flags[prot->num_atoms];
    make_asa_list(prot,-1,-1,-1,flags);
    
    double probe = 1.4;
    int ndiv = 3;
    double ext_radius = 1.4;
    
    CntmultList cnt_list[prot->num_res];
    
    double score = asa_score(prot,flags,probe,ndiv,wat_list,wlist,ext_radius);
    printf("Score %lf",score);
    
    */
     //decorate_box(argv[1],argv[2],argv[3]);

    
    /*
    decorate_box(prot,
    177.560, -177.347,  183.661,   
    -170.185,-171.188,  181.350,  
    -172.754, -171.844, 177.706,
    -169.012, -167.163,  181.003
    );
    */
    /*
    LinusProtein *prot = malloc(sizeof(LinusProtein));
   
    protein_from_pdb(argv[1],prot);
    Atom3d *atoms = prot->atoms;
  
    FILE *file = fopen(argv[2],"r");
    char buff[BUFSIZ];
    int k=0;
    while (fgets(buff, sizeof buff, file) != NULL)
    {
        sscanf(buff, "%lf %lf %lf\n",&atoms[k].x,&atoms[k].y,&atoms[k].z);
        k+=1;
    }
    int m = bumpcheck(atoms, prot->res_fai, 1, prot->num_res-2, 2, prot->num_res-1);
    if (m) 
    {
        printf("Bump %d\n",m);
    }
    else
    {
        printf("No Bump %d\n",m);        
    }     
*/
  
  

/*  
    HbList *hblist;
    int i,j,size=0,max_acp_size=0;
       
    int num_hb = make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&size,&max_acp_size);
    hblist=malloc(size*sizeof(HbList));
   
    for(i=0;i<size;i++)
    {
        hblist[i].acps = malloc(max_acp_size*sizeof(Acps));
    }

    make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&num_hb,&max_acp_size);
    for(i=0;i<num_hb;i++)
    {
        for(j=0;j<hblist[i].acps_size;j++)
        {
            printf("Hblist %s %s %s %s %s\n",prot->atoms[hblist[i].donor].name,prot->atoms[hblist[i].da1].name,prot->atoms[hblist[i].da2].name,
            prot->atoms[hblist[i].acps[j].acp].name, prot->atoms[hblist[i].acps[j].acp1].name);
        }
    }
    
    
    //Linus_sim sim;
    
    //create_linus_sim(&sim,argv[1]);
    
    //run_sim(&sim);
    
    //double w[18] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0};
    //update_weights(&sim,w);
    
    //read_weights_from_file(&sim,"wt.txt",4);
    
    //char *type="helix"; 
    
    //normalize_weights(&sim);
    
    
    /*
    int i;
   
    for (i=0;i<sim.protein.num_res;i++)
    {
        printf("Main %lf\n",sim.protein.res_helix_wt[i]);
    }
    
    /*
    LinusFile lf;
    create_linusfile(&lf, argv[1]);
    FILE *f;
    get_logfile(&lf,f);
    printf("Trying to write to file\n");
    fprintf(f,"Writing to log file works!\n");
    */
    
    
    //double s[6];
    
    //gen_seeds(s);
    
    //smrt_pmeso(prot,"pmesowt.txt");
    //double b[3]={1.0,2.0,3.0};
    //double x = _distance(&prot->atoms[0],b);
    //printf("_distance %lf\n",x);
    
    //double rmsd;
    //int x = local_move(prot,0,3,s,&rmsd);
    //printf("Success %d RMSD %lf\n",x,rmsd);
    /*
    int i;
    for(i=0;i<=6;i++)
    {
        printf("%lf ",s[i]);
    }
    printf("\n");
    
    int x = choice (6,s);
    printf("Choice output %d \n",x);
    
    double m = uniform(1.0,20.0,s);
        printf("Uniform - %lf \n",m);
    */
    //
    
    /*
    char *sst[19];
    int i;
    for (i = 0; i < 19; i++)
    {
        sst[i] = (char *)malloc(sizeof(char)*2);
    }
    
    
    rc_ss(prot,sst);
    */
        
    
    /*
    char **codes;
    int i; 
    //printf("Start allocating\n");
    codes = (char **)malloc(prot->num_res * sizeof(char *));
    
    for (i = 0; i < prot->num_res; i++)
    {
        codes[i] = (char *)malloc(sizeof(char)*3);
    }
    //printf("Finishde allocating\n");
    
    rc_codes(prot,codes);
   

    for (i=0;i<prot->num_res;i++)
    {
        printf("Code for res %s is %s\n",prot->res_names[i],codes[i]);
    }
    
    //get_rc_code(-180,-165,code);
    //res_rc(150.0,165.0,0.0,code);
    */
    
    
    /*
    HbList *hblist;
    int i,j,size=0,max_acp_size=0;
    
   
    make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&size,&max_acp_size);
   
    hblist=malloc(size*sizeof(HbList));
   
    for(i=0;i<size;i++)
    {
        hblist[i].acps = malloc(max_acp_size*sizeof(Acps));
    }

    make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&size,&max_acp_size);
    int num = calc_hbs_score_local(prot,hblist,size,0);
    printf("Num %d\n",num);
    */
    
    
    
    /*
    int numint_loos, numvirt_loos, bbtot,numbb;
    CntmultList solv_list[prot->num_res];
    double *wat_list;
    
    FILE *file;
    
    int wlist = mk_virt_bb_cntmult_loosHbd(prot,hblist,size,0,0,file,0,1.25,&numint_loos,&numvirt_loos,&bbtot,&numbb,wat_list,0,&solv_list);
    wat_list = (double *) malloc(wlist*sizeof(double));
    numint_loos= numvirt_loos=bbtot=numbb=0;
    mk_virt_bb_cntmult_loosHbd(prot,hblist,size,0,0,file,0,1.25,&numint_loos,&numvirt_loos,&bbtot,&numbb,wat_list,wlist,&solv_list);
    

    
    int flags[prot->num_atoms];
    make_asa_list(prot,-1,-1,-1,flags);
    
    double probe = 1.4;
    int ndiv = 3;
    double ext_radius = 1.4;
    
    CntmultList cnt_list[prot->num_res];
    
    double score = print_asa_score(prot,flags,probe,ndiv,wat_list,wlist,ext_radius,&cnt_list ,"asa_out");
    printf("Score %lf",score);
    
    */
   
    
    
    
   /*
   double score = calc_chasa_solv(prot,hblist,size);
   printf("Score %lf\n",score);
 
    */
    
    
    
    
    
    
    
    
    
    
    
    
    //void print_hbs_casa(LinusProtein *p, int * flags, double probe, int ndiv,char *filename)

    /*
    for (i=0;i<=prot->num_res;i++)
    {
        printf("i %d fai %d\n",i, prot->res_fai[i]);
    }
    */
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
    Atom3d *satoms = malloc(prot->num_atoms * sizeof(Atom3d));
    int size = set_solvation_parameters(prot, 1.4, satoms);
    
    
    
    double asa = geocav(size,satoms,4);
    printf("asa %lf\n",asa);
    
    double sasa = solvation_score(prot);
    
    printf("Sasa %lf\n",sasa);
    
    
    
    int flags[prot->num_atoms];
    make_asa_list(prot,-1,-1,-1,flags);
    print_hbs_casa(prot,flags,-1.0,-1,"hbs_casa.txt");
    
    
    /*
  
 

    
    /*
    AtomList *hbtab_loos;
    int hbtab_size = make_hbtab_loos(prot,hblist,size,0,hbtab_loos);
    
    hbtab_loos=malloc(hbtab_size*sizeof(AtomList));
    
    make_hbtab_loos(prot,hblist,size,hbtab_size,hbtab_loos);
    
    for(i=0;i<hbtab_size;i++)
    {
        printf("Hbtab %s %d %s %d  %lf\n",hbtab_loos[i].a1.name,hbtab_loos[i].a1.resnum,hbtab_loos[i].a2.name,hbtab_loos[i].a2.resnum,hbtab_loos[i].value);
    }
    
    int numint_loos=0;
    for(i=0;i<prot->num_atoms;i++)
    {
         int x = _ishbonded_loos(prot,&prot->atoms[i],hbtab_loos,hbtab_size,&numint_loos);
         printf("%s %d %d %d\n",prot->atoms[i].name,prot->atoms[i].resnum,x,numint_loos);
    }
    
        // mk_virt_bb_cntmult_loosHbd
    int i;
    HbList *hblist;
    int j,size=0,max_acp_size=0;
    
   
    make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&size,&max_acp_size);
   
    hblist=malloc(size*sizeof(HbList));
   
    for(i=0;i<size;i++)
    {
        hblist[i].acps = malloc(max_acp_size*sizeof(Acps));
    }

    make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&size,&max_acp_size);
   
 
    FILE *file = fopen("water.pdb","w+");
    
    int numint,numvirt,bbtot,numbb;
    double *wat_list;
    int ctr_wat_list =0;
    CntmultList cntmult_list[prot->num_res];
    
    int wlist = mk_virt_bb_cntmult_loosHbd(prot,hblist,size,0,0,file,1,1.25,&numint,&numvirt,&bbtot,&numbb,wat_list,0,&cntmult_list);
    
    printf("\n\nGot the size of wlist %d. Now populating water list\n",wlist);
    
    wat_list = (double *) malloc(wlist*sizeof(double));
    mk_virt_bb_cntmult_loosHbd(prot,hblist,size,0,0,file,1,1.25,&numint,&numvirt,&bbtot,&numbb,wat_list,wlist,&cntmult_list);
    
    
    
    
    // mk_virt_bb_cntmult_loosHbd end
    
    
    
    
    
    
    
    /*
    int flags[prot->num_atoms];
    make_asa_list(prot,-1,-1,-1,flags);
    
    int i;
     for (i=0;i<prot->num_atoms;i++)
    {
        printf("Flag for atom %s %d is %d\n",prot->atoms[i].name,prot->atoms[i].resnum,flags[i]);
    }
    
    double *extra_coords;
    double score = print_asa_score(prot, flags, 1.4,3,extra_coords,1.4,"asa_out.txt");
    printf("score %lf\n",score);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
    int num = add_waters(prot);
    
    int i;
    for(i=0;i<num;i++)
    {
        printf("%s %s %s %lf %lf %lf\n", prot->waters[i].first_parent->name,prot->waters[i].second_parent->name,prot->waters[i].third_parent->name,prot->waters[i].fp_distance,prot->waters[i].sp_angle,prot->waters[i].tp_torsion);
    }
    
    
    
    
    
    
    /*
    char **codes;
    int i; 
    //printf("Start allocating\n");
    codes = (char **)malloc(prot->num_res * sizeof(char *));
    
    for (i = 0; i < prot->num_res; i++)
    {
        codes[i] = (char *)malloc(sizeof(char)*3);
    }
    //printf("Finishde allocating\n");
    rc_codes(prot,codes);
   

    for (i=0;i<prot->num_res;i++)
    {
        printf("Code for res %d is %s\n",prot->res_names[i],codes[i]);
    }
    
    //get_rc_code(-180,-165,code);
    //res_rc(150.0,165.0,0.0,code);
    /*
    int i;
    AtomList *vlist;
    int size = make_coul_list(prot,vlist,0);
    vlist=malloc(size*sizeof(AtomList));
    make_coul_list(prot,vlist,size);
    
   // for(i=0;i<size;i++)
    //{
      //  printf("A1 %s A2 %s Sigma %lf\n",vlist[i].a1.name,vlist[i].a2.name,vlist[i].value);
    
   // }
    
    
    double score = coul_score(prot,vlist,size,1.0);
    printf("score %lf\n",score);
    
    
    
    
    /*
    
    
    ----------------hbond list test--------------------------------------------------------------
    HbList *hblist;
    int i,j,size=0,max_acp_size=0;
    
    printf("Getting size\n");
    make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&size,&max_acp_size);
   
    printf("Size %d Max_acp_size %d\n",size,max_acp_size); 
    hblist=malloc(size*sizeof(HbList));
   
    
    for(i=0;i<size;i++)
    {
        hblist[i].acps = malloc(max_acp_size*sizeof(Acps));
    }
    
    printf("\n\nSize of hblist %d\n",size);
    
    printf("\n\nPopulating List\n\n");
    make_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,0.5,1.5,3.0,130.0,0.5,1,1,&size,&max_acp_size);
     
    print_hbond_list(prot,hblist,size,"hbondlist.txt"); 


    printf("\nReading list\n");
    for(i=0;i<size;i++)
    {
        printf("%s %d %s %s\n",hblist[i].donor.name,hblist[i].donor.resnum, hblist[i].da1.name, hblist[i].da2.name);
        for (j=0;j<hblist[i].acps_size;j++)
       {
            //printf("Size of acp is %d %d %d\n",hblist[i].acps_size,i,j);    
            printf("    %s %s \n",hblist[i].acps[j].acp.name,hblist[i].acps[j].acp1->name);
       }
    }
    --------------End hbond list test--------------------------------------------------------------
    ----------------Chasa test--------------------------------------------------------------
    HbList *hblist;
    int i,j,size=0,max_acp_size=0;
    
    printf("Getting size\n");
    make_chasa_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,3.0,130.0,&size,&max_acp_size);
   
    printf("Size %d Max_acp_size %d\n",size,max_acp_size); 
    hblist=malloc(size*sizeof(HbList));
   
    
    for(i=0;i<size;i++)
    {
        hblist[i].acps = malloc(max_acp_size*sizeof(Acps));
    }
    
    printf("\n\nSize of hblist %d\n",size);
    
    printf("\n\nPopulating List\n\n");
    make_chasa_hbond_list(prot,hblist,2,6,3.0,1.5,130.0,3.0,130.0,&size,&max_acp_size);
     
    

    printf("\nReading list\n");
    for(i=0;i<size;i++)
    {
        printf("%s %d %s %s\n",hblist[i].donor.name,hblist[i].donor.resnum, hblist[i].da1.name, hblist[i].da2.name);
        for (j=0;j<hblist[i].acps_size;j++)
       {
            //printf("Size of acp is %d %d %d\n",hblist[i].acps_size,i,j);    
            printf("    %s %s \n",hblist[i].acps[j].acp.name,hblist[i].acps[j].acp1->name);
       }
    }
    
    ----------------End Chasa test--------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    AtomList *clist;
    int size = make_contact_list(prot,0,4,1.0,clist,0);
    clist=malloc(size*sizeof(AtomList));
    make_contact_list(prot,0,4,1.0,clist,size);
    int x;
    for(x=0;x<size;x++)
    {
        printf("%s %d %s %d %lf %d\n",clist[x].a1.name, clist[x].a1.resnum, clist[x].a2.name, clist[x].a2.resnum,clist[x].value, x);    
    }
    
    
    double cscore = contact_score(clist,size,1.0);
    printf("Contact score %lf\n", cscore);
    cscore = print_contact_score(prot,clist,size,1.0,"contact_scores.txt");
    printf("Total contact score after printing %lf\n", cscore);
    
    
    
    
    double score;
    print_compaction_score(prot,1.0,"compact_score.txt",&score);
    printf("Compaction score %lf",score);
  
  
    double c1[3]={1.0,2.0,3.0};
    double c2[3]={1.0,2.0,3.0};
    int x = _close(c1,c2,2.0);
    printf("Close %d",x);
    
    
    double coordsfp[3]={1.0,2.0,3.0};
    double coordssp[3]={4.0,5.0,6.0};
    double coordstp[3]={-7.0,8.0,-9.0};
    int zmindx[3]={1,2,3};
    double zmvalu[3]={1.1,2.2,3.3};
    double atcoor[3];
    _ztox(coordsfp,coordssp,coordstp,zmvalu,zmindx,atcoor);
    printf("%lf %lf %lf\n", atcoor[0],atcoor[1],atcoor[2]);
    
    
     double x_norm,y_norm,z_norm;
    
    _norm(-5.0,4.0,-3.0,1.0,2.0,3.0,&x_norm,&y_norm,&z_norm);
    
    int flags[prot->num_atoms];
    make_asa_list(prot,-1,-1,-1,flags);
    
    int i;
     for (i=0;i<prot->num_atoms;i++)
    {
        printf("Flag for atom %s %d is %d\n",prot->atoms[i].name,prot->atoms[i].resnum,flags[i]);
    }
    
    /*
    AtomList *clist;
    int size = make_coul_list(prot,clist,0);
    clist=malloc(size*sizeof(AtomList));
    make_coul_list(prot,clist,size);
    
    int x;
    for(x=0;x<size;x++)
    {
        printf("%s %d %s %d %lf %d\n",clist[x].a1.name, clist[x].a1.resnum, clist[x].a2.name, clist[x].a2.resnum,clist[x].value, x);    
    }
    
    
        int i;
    for(i=0;i < num_vlist;i++)
    {
          printf("%s %d %s %d %lf\n",vlist[i].a1.name,vlist[i].a1.resnum,vlist[i].a2.name,vlist[i].a2.resnum,vlist[i].value);
    }
    
    
    
    
   
        //double coord[3];
         //printf("%s %d ", prot->atoms[i].name, prot->atoms[i].resnum);
        //get_coords(&prot->atoms[i], coord);
        //printf("X %lf Y %lf Z %lf\n",coord[0],coord[1],coord[2]);
        
       /* 


        int i;
        for(i=0;i < prot->num_atoms;i++)
        {
              
        }
               
        double sqdist1;
        double dist1;
        distance(&prot->atoms[i],&prot->atoms[i+1],&dist1);
        sqdistance(&prot->atoms[i],&prot->atoms[i+1],&sqdist1);
        printf("%s %d %lf %lf %lf %lf\n",prot->atoms[i].name, prot->atoms[i].resnum, dist1, prot->atoms[i].fp_distance, prot->atoms[i].sp_angle, prot->atoms[i].tp_torsion);
        
        double dx =   prot->atoms[i].x - prot->atoms[0].x;   
        double dy =   prot->atoms[i].y - prot->atoms[0].y;
        double dz =   prot->atoms[i].z - prot->atoms[0].z;
        int y =   close(&prot->atoms[0], &prot->atoms[i], -1.0);
        printf(" %d %lf %lf %lf %lf", y, dx*dx, dy*dy, dz*dz, dx*dx+dy*dy+dz*dz);
        if(y)
            printf("Too close!!!\n");
        else
            printf("Far apart!\n");
            
        prot->atoms[i].fp_distance_tmp = 10.0;
        prot->atoms[i].sp_angle_tmp = 10.0;
        prot->atoms[i].tp_torsion_tmp = 10.0;
        
        restore_coords(&prot->atoms[i]);
        ztox(prot->atoms, 0,prot->num_atoms);
        double rg1;
        radius_of_gyration(prot->atoms,0,prot->num_atoms,&rg1);
        rotatex(&prot->atoms[1],30.0);
    
        int is_bump = bump_check_intra(prot->atoms,prot->res_fai,i);
        printf("%s %d \n",prot->res_names[i],is_bump);     
         print_all_bumps(prot,"test_bumps.txt");
         
        char *atom_name = " N  ";
        int is_H = is_hydrophobic(atom_name);
        printf("Is hydrogen %d\n", is_H);
        
        int minres, maxres;
    get_res_extents(prot,&minres,&maxres);
    printf("%d %d\n",minres,maxres);
     Atom3d *scatoms;
    int num_sc = get_scacceptor(prot->atoms,0,8,"GLU",scatoms);
    
   
    int i;
    for(i=0;i < num_sc;i++)
    {
           printf("Acceptor atoms %s\n",scatoms[i].name);
    }
      

            Atom3d *satoms;
            set_solvation_parameters(prot,1.4,satoms);

      
    } 
    
    
    AtomList * vlist;
    int vlist_size = make_vdw_list(prot,vlist,0);
    vlist=malloc(vlist_size*sizeof(AtomList));
    int num_vlist = make_vdw_list(prot,vlist,vlist_size);
    
    int x;
    for(x=0;x<num_vlist;x++)
    {
        printf("%s %d %s %d %lf %d\n",vlist[x].a1.name, vlist[x].a1.resnum, vlist[x].a2.name, vlist[x].a2.resnum,vlist[x].value, x);    
    }
    
    */
        

    
    return 0;
}   



