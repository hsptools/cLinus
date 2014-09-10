typedef struct
{
    double beta;
    int numcycle;
    int numsave;
    int trials;
    int updateweights;
    int updatepmeso;
    int reset_conformation;
    int reset_emin_struc;
    int use_lock_ss;
    int use_msd_focus;
    int use_msd_sampling;
    int use_ms_sampling;
    int use_pms_sampling;
    int use_local_mov;
    int use_dist_constraints;
    int write_rc;

}Simpar;


typedef struct
{
    int use_hbond; //: 0,
    int use_sidechain_hbond; // 0,
    double hbond_distance; // 3.0,
    double hbond_probe; // 1.5,
    double hbond_score_short; // 0.5,
    double hbond_score_long; // 1.5,
    double hbond_torsion; // 130.0,
    double sidechain_hbond_distance; // 3.0,
    double sidechain_hbond_score; // 0.5,
    double sidechain_hbond_torsion; // 130.0,
    int hbond_winmin; // 2,
    int hbond_winmax; // 6,
    int use_hbs_local; // 0

}Hbdpar;

typedef struct
{
    int use_contact; // 0,
    double contact_probe; // 2.8,
    int contact_winmin; // 2,
    int contact_winmax; // 6,
    double contact_scale; // 1.0,
    int use_softsph; // 0,
    int use_LJ; // 0
 
}Conpar; 

typedef struct
{
    int use_torsion; // 0,
    double torsion_penalty; // 1.0,
}Torpar;

typedef struct
{
    int use_compact; // 0,
    double compact_expect; // 1.0
}Compar;

typedef struct
{
    int use_chasa; // 0
}Chasa;

typedef struct
{
    int use_solvation; // 0,
    double solvation_probe; // 1.4,
    double cavity_term; // 0.0,
    int solvation_winmax; // 6
}Slvpar;






















