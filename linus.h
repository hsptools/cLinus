#include "Atom3d.h"
#include "LinusProtein.h"
#include "AtomList.h"
#include "HbList.h"
#include "DcList.h"
#include "CntmultList.h"
#include "LinusFile.h"
#include "parameters.h"
#include "Linus_sim.h"



void find_native_like( double nat_phi1, double nat_psi1, double nat_ome1,
                       double nat_phi2, double nat_psi2, double nat_ome2,
                       double nat_phi3, double nat_psi3, double nat_ome3,
                       double nat_phi4, double nat_psi4, double nat_ome4, 
                       double rangep,char *allowed, char *prefix);

double set_sy(double psi);

double linus_sy(double psi);
double out_sy(double psi);
double _tor(Atom3d *a1, Atom3d *a2, Atom3d *a3, Atom3d *a4);
void print_conf(char *fmt, LinusProtein *prot, FILE *ofile,int num_steric_hbs);
double res_rmsd(Atom3d *a,Atom3d *b,Atom3d *c,Atom3d *d,double *e,double *f,double *g,double *h);
double _distance(Atom3d *a,double *b);

double random1(double *s);
int choice(int size, double *s);
double uniform(double start, double end, double *s);

int mk_virt_bb_cntmult_loosHbd(LinusProtein *p,HbList *hblist,int hblist_size, int res, int tryall,FILE *watpdb,int watpdb_write,double ext_rad,
int *numint_loos, int *numvirt_loos, int *bbtot, int *numbb, double *wat_list, int wat_list_size, CntmultList *cntmult_list);

double asa_score(LinusProtein *p, int  *flags, double probe, int ndiv, double *ext_coords, int numext, double ext_radius);
double asa_eval(Atom3d * atoms, int numatm, double * data, double *coords, int *flags,double probe, double ext_rad, int use_data, int numext, int ndiv);

int  find_neighbors(Atom3d *atoms, double *coords,int *flags,double probe, 
double ext_rad,int k, int numatm, int numext, int *nghlst);
int tri_buried(Atom3d* atoms, double *coords, int numext, double tx, double ty, double tz, 
double probe, double ext_rad, int *nghlst, int numngh,int *nghstrt);
double print_asa_score(LinusProtein *p,int  *flags, double probe, int ndiv, double *ext_coords, int numext, double ext_radius, CntmultList *solv_list,char *filename);

double calc_chasa_solv2(LinusProtein *p, int numint_loos, int numvirt_loos, int btot, int numbb, double *wat_list, int wlist,CntmultList *solv_list);

double solvation_score(LinusProtein * p);

double geocav(long ncor, Atom3d *atoms, int winmax);

double distance(Atom3d *a1, Atom3d *a2);

double sqdistance(Atom3d *a1, Atom3d *a2);

double torsion(Atom3d *a1, Atom3d *a2, Atom3d *a3, Atom3d *a4);

int close(Atom3d *a1, Atom3d *a2, double d);

double angle(Atom3d *a1, Atom3d *a2, Atom3d *a3);

void _norm(double x1, double y1, double z1, double x2, double y2, double z2, double * x_norm, double *y_norm, double *z_norm);


void get_coords(Atom3d *self, double * coords);

void get_atom_resid(LinusProtein *prot, int resid, char *aname, Atom3d *atom);

void get_res_extents(LinusProtein *p, int *minres, int *maxres);

double get_sasa_rad(char *rname, char *aname);


double calc_chasa_solv(LinusProtein *p, HbList *hblist, int hblist_size);
int make_contact_list(LinusProtein * p, int wmin, int wmax, double cscale, AtomList * clist, int size);

int make_coul_list(LinusProtein * p, AtomList * qlist, int size);
double coul_score(LinusProtein *p,AtomList *qlist, int size, double scale);

int make_vdw_list(LinusProtein * p, AtomList * vlist, int size);
double soft_LJ_score(LinusProtein *p, AtomList *vlist, int size, double scale);
double LJ_score(LinusProtein *p, AtomList *vlist, int size, double scale);
double soft_sph_score(LinusProtein *p, AtomList *vlist, int size);
double soft_debump_score(LinusProtein *p, AtomList *vlist, int size);
double print_soft_sph_score(LinusProtein *p, AtomList *vlist, int size, char*filename);
double print_LJ_score(LinusProtein *p, AtomList *vlist, int size, double scale,char *filename);

double dist_const_score(LinusProtein *p,DcList *dclist, int size);

int make_asa_list(LinusProtein *p, int hphob, int hphil, int hydrogens, int *flags);


double radius_of_gyration(Atom3d *atoms, int start, int end);
void geocenter(Atom3d *atoms,int start, int end, double *x, double *y, double *z );

int bumpcheck(Atom3d *atoms, int *fai, int ires1, int ires2, int jres1, int jres2);
void print_all_bumps(LinusProtein *prot, char *filename);


void get_atom_with_name(Atom3d *atom, Atom3d *atoms, int start, int end, char *aname);

double contact_score(LinusProtein *p,AtomList * clist,int size, double cprobe);
double print_contact_score(LinusProtein *p, AtomList *clist, int size,double cprobe, char *filename);

double hbond_score (LinusProtein *p,HbList *hblist, int nhb);

double torsion_score(LinusProtein *p, double tp);



double compaction_score(LinusProtein *p, double ce);
void print_compaction_score(LinusProtein *p, double cp, char *filename, double *score);

int make_hbond_list(LinusProtein * p, HbList * hblist, int wmin, int wmax,double hbdist, 
                    double hbprobe, double hbtor, double hbenes,
                    double hbenel, double shbdist,  double shbtor, double shbene, 
                    int use_hbond, int use_schbondn, int *size,int *max_acp_size);

int make_chasa_hbond_list(LinusProtein * p, HbList * hblist, 
                     int wmin, int wmax, 
                     double hbdist, double hbprobe, 
                     double hbtor,  double shbdist, 
                     double shbtor, int *size, int *max_acp_size);                    
                    
                    
/*

default_hbdpar = {
    'use_hbond': 0,
    'use_sidechain_hbond': 0,
    'hbond_distance': 3.0,
    'hbond_probe': 1.5,
    'hbond_score_short': 0.5,
    'hbond_score_long': 1.5,
    'hbond_torsion': 130.0,
    'sidechain_hbond_distance': 3.0,
    'sidechain_hbond_score': 0.5,
    'sidechain_hbond_torsion': 130.0,
    'hbond_winmin': 2,
    'hbond_winmax': 6,
    'use_hbs_local' : 0
}

*/


