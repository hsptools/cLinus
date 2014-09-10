
/*
"""Class to describe a LINUS simulation object

    Instantiation

        An instance is created by calling the class constructor with
        two arguments - name of the pdb file and a random number generator
        The pdb file must reside in a directory to which the USER has
        write priviledges.
       

    Attributes

        o *protein* - instance of LinusMol - the protein that is being
        simulated

        o *simpar* - dictionary - simulation parameters. The various
        parameters are

            i.    *beta* - float-  reciprocal of the simulation temperature
            
            ii.   *numcycle* - integer - number of cycles of simulation to run
            
            iii.  *numsave*  - integer - frequency with which structures are
                  saved to disk
                              
            iv.   *reset_conformation* - boolean - if true at the end of
                  the simulation the molecules' conformation is reset
                  to the one defined in the input pdb file
                  
            v.    *update_weights* - boolean - if true at the end of the
                  simulation the sampling weights of the residues are
                  updated to the observed distribution in the various
                  secondary structural states during the course of the
                  simulation
                  
                  
                  
        o *hbdpar* - hydrogen bonding energy parameters.

            i.   *use_hbond* - boolean - if true backbone - backbone
                 hydrogen bonds are enabled

            ii.  *use_sidechain_hbond* - boolean - if true sidechain
                 acceptor to backbone donor hydrogen bonds are
                 enabled

            iii. *hbond_distance* - float - optimal distance between
                 backbone donor and backbone acceptor

            iv.  *hbdond_torsion* - float - minium torsion angle between
                 the backbone acceptor and the backbone donor and it's two
                 antecedent atoms( 'C' of previous residue and 'CA' of the
                 residue

            v.   *hbond_probe* - float - distance over which the hydrogen
                 bond energy for backbone donor to backbone acceptor
                 bonds will be scaled to zero from the maximum value

            vi.  *hbond_score_short* - float - backbone donor to backbone
                 acceptor hydrogen bond energy for residues separations of
                 5 or less

            vii. *hbond_score_long* - float - backbone donor to backbone
                 acceptor hydrogen bond energy for residues separations of
                 6 or less

            viii. *sidechain_hbond_distance* - float - optimal distance between
                  backbone donor and sidechain acceptor

            ix.  *sidechain_hbond_torsion* - float - minium torsion angle
                 between the sidechain acceptor and the backbone donor and
                 it's two antecedent atoms ( 'C' of previous residue and 'CA' 
                 of the residue)

            x.   *sidechain_hbond_score* - float - energy of a sidechain
                 to backbone hydrogen bond

            xi.  *hbond_winmin* - integer - minimum separation between
                 two residues that can participate in backbone to backbone
                 hydrogen bonds

            xii. *hbond_winmax* - integer - maximum  separation between
                 two residues that can participate in backbone/sidechain
                 acceptor to backbone donor hydrogen bonds

        o *conpar* - contact energy parameters

            i.   *use_contact* - boolean - if true contact energy
                 calculations are enabled

            ii.  *contact_winmin* - integer - minimum sequence separation
                 between pairs of residues that participate in contact
                 energy calculations

            iii. *contact_winmax* - integer - maximum sequence separation
                 between pairs of residues that participate in contact
                 energy calculations

            iv.  *contact_probe* - float - distance over which the
                 contact energy is scaled to zero from the maximum value

            v.   *contact_scale* - float - scaling factor to apply to
                 the contact energy

        o *torpar* - torsion penalty score parameters

            i.   *use_torsion* - boolean - if true torsion energy calculation
                 is enabled

            ii.  *torsion_penalty* - float - penalty for adoption positive
                 phi values - for GLY this is a reward instead of a penalty

        o *files* - instance of LinusFiles - class to create various files
          for logging etc. during a simulation

        o *generator* - object - a random number generator

        o *registered functions* - list - list of functions to be run
          before, during and after a simulation
          

*/
typedef struct
{
    char *func;
    int when;
}Register;


typedef struct
{
    char *pdbfile;
    LinusProtein protein;
    Simpar simpar;
    Hbdpar hbdpar;
    Conpar conpar;
    Compar compar;
    Torpar torpar;
    Slvpar slvpar;
    Chasa chasa;
    int *movelist;
    LinusFile files;
    double s[6];
    HbList *hblist;
    int num_hb;
    AtomList *vlist;
    int num_v;
    AtomList *conlist;
    int num_con;
    DcList *dclist;
    int num_dc;
    Register registered_functions[100];
    int max_regfun;
    int trial;
    int trials;
    double *swt; // list of lists - each list item is a list of length 6 containing, respectively, the helix, strand, turn1, turn2,coil and PII weights.
    

}Linus_sim;


