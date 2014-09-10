#include <stdio.h>
#include "linus.h"

void default_simpar(Simpar *simpar)
{
    simpar->beta= 2.0;
    simpar->numcycle= 1000;
    simpar->numsave= 10;
    simpar->trials= 1;
    simpar->updateweights= 0;
    simpar->updatepmeso= 0;
    simpar->reset_conformation= 0;
    simpar->reset_emin_struc= 0;
    simpar->use_lock_ss= 0;
    simpar->use_msd_focus= 0;
    simpar->use_msd_sampling= 0;
    simpar->use_ms_sampling= 0;
    simpar->use_pms_sampling= 0;
    simpar->use_local_mov= 0;
    simpar->use_dist_constraints= 0;
    simpar->write_rc= 0;

}

void default_hbdpar(Hbdpar *hbdpar)
{
    hbdpar->use_hbond= 0;
    hbdpar->use_sidechain_hbond= 0;
    hbdpar->hbond_distance= 3.0;
    hbdpar->hbond_probe= 1.5;
    hbdpar->hbond_score_short= 0.5;
    hbdpar->hbond_score_long= 1.5;
    hbdpar->hbond_torsion= 130.0;
    hbdpar->sidechain_hbond_distance= 3.0;
    hbdpar->sidechain_hbond_score= 0.5;
    hbdpar->sidechain_hbond_torsion= 130.0;
    hbdpar->hbond_winmin= 2;
    hbdpar->hbond_winmax= 6;
    hbdpar->use_hbs_local = 0;

}

void default_conpar(Conpar *conpar)
{
    conpar->use_contact= 0;
    conpar->contact_probe= 2.8;
    conpar->contact_winmin= 2;
    conpar->contact_winmax= 6;
    conpar->contact_scale= 1.0;
    conpar->use_softsph= 0;
    conpar->use_LJ= 0;
}


void default_torpar(Torpar *torpar)
{
  torpar->use_torsion= 0;
    torpar->torsion_penalty= 1.0;

}


void default_compar(Compar *compar)
{
 compar->use_compact= 0;
    compar->compact_expect= 1.0;

}


void default_chasa(Chasa *chasa)
{
    chasa->use_chasa= 0;

}

void default_slvpar(Slvpar *slvpar)
{

  slvpar->use_solvation= 0;
    slvpar->solvation_probe= 1.4;
    slvpar->cavity_term= 0.0;
    slvpar->solvation_winmax= 6;
}

// *set_simulation_parameters* - set some or all of the simulation parameters



void create_linus_sim(Linus_sim *sim, char *pdbfile)
{
    protein_from_pdb(pdbfile, &sim->protein); 
    default_simpar(&sim->simpar);
    default_hbdpar(&sim->hbdpar);
    default_conpar(&sim->conpar);
    default_compar(&sim->compar);
    default_torpar(&sim->torpar);
    default_slvpar(&sim->slvpar);
    default_chasa(&sim->chasa);
    sim->trial =0;
    sim->trials =0;
    sim->max_regfun =0;
    gen_seeds(sim->s);
    create_linusfile(&sim->files,pdbfile);
    //self.hblist = self.conlist = self.dscnstrlist = None
    //self.vlist = None
    //self.generator = generator
    //self.registered_functions = []
}

