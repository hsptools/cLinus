typedef struct
{
// Attributes
char name[5];
char filename[30];                 //string - name of file from which the molecule information was read
Atom3d *atoms;                 //list of atoms in the molecule. each item in this list is an instance of *Atom3d*
int *satoms;
int num_atoms;                  //number of atoms in molecule
int num_res;                     //number of residues in molecule
char **res_names;                 //list - name of each residue in molecule
int *res_num;                   //list - number of each residue in input PDB file 
int *res_fai;                   //list - index of the the first atom of each residue (this is horrendously named)
int *phi_atoms;              //list - contains a reference to the 'C' atom of each residue
int *psi_atoms;               //list - contains a reference to the 'O' atom of each residue
int *omega_atoms;             //list - contains a reference to the 'CA' atom of each residue
int *chi_atoms;
int num_chi_atoms;
Atom3d *chi_torsions;           //list of lists - references to the atoms corresponding to the various chi torsions of each residue
Atom3d *waters;
double *res_helix_wt;    //list of floats - helix sampling weights for each residue
double *res_strand_wt;   //list of floats - strand sampling weights for each residue
double *res_coil_wt;     //list of floats - coil sampling weights for each residue
double *res_turn1_wt;    //list of floats - turn1 sampling weights for each residue
double *res_turn2_wt;    //list of floats - turn2 sampling weights for each residue
double *res_PII_wt;      // list of floats
double *orig_fp_distance;
double *orig_sp_angle;
double *orig_tp_torsion;
double *fp,*sp,*tp;
double *res_f_mean;
double *res_f_sd;
double *res_y_mean;
double *res_y_sd;
double *res_acprat;
double pmeso_grid[144];
double **pmeso_wt;
char **mesostates;                //list of strings - allowed mesostates for each residue. *Note* - this attribute is present only if mesostates have been specified
}LinusProtein;

/*
 """Class to describe a protein for simulation in LINUS.

    Instantiation

       Instantiation of this class requires one argument - the
       name of the file containing a description of the molecule
       in PDB format.

    Class Attributes

       o *filename* - string - name of file from which the molecule
       information was read

       o *atoms* - list - list of atoms in the molecule. each item
       in this list is an instance of *Atom3d*

       o *num_atoms* - integer - number of atoms in molecule

       o *num_residues* - integer - number of residues in molecule

       o *residue_names* - list - name of each residue in molecule

       o *res_pdb_number* - list - number of each residue in input PDB file 

       o *residue_first_atom_indices* - list - index of the the
       first atom of each residue (this is horrendously named)

       o *phi_atoms* - list - contains a reference to the 'C' atom
       of each residue

       o *psi_atoms* - list - contains a reference to the 'O' atom
       of each residue

       o *omega_atoms* - list - contains a reference to the 'CA' atom
       of each residue

       o *chi_torsions* - list of lists - references to the atoms
       corresponding to the various chi torsions of each residue

       o *residue_helix_weights* - list of floats - helix sampling weights
       for each residue

       o *residue_strand_weights* - list of floats - strand sampling weights
       for each residue

       o *residue_coil_weights* - list of floats - coil sampling weights
       for each residue

       o *residue_turn1_weights* - list of floats - turn1 sampling weights
       for each residue

       o *residue_turn2_weights* - list of floats - turn2 sampling weights
       for each residue

       o *mesostates* - list of strings - allowed mesostates for each
       residue. *Note* - this attribute is present only if mesostates
       have been specified


*/


