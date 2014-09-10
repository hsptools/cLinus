/*
 o *GLY_COIL* - the possible coil values for glycine residues.
     This is a tuple of tuples, where each tuple lists the
     range of phi, psi and omega respectively

   o *PRO_COIL* - possible coil values for proline

   o *NPG_COIL* - possible coil values for all residues that are not
     proline or glycine.

   o *PRO_PRO_TURN* - possible turn values for two successive proline
     residues.  This is a tuple of tuples. Each individual tuple is of
     length 2, with the first item being a tuple of length 3 listing
     the range of the phi, psi and omega values for the residue at
     position i+1 of a beta-turn and the second item listing the
     values for the i+2 residue.

   o *PRO_X_TURN* - possible turn values for  two successive residues
     in which the first residue is a proline

   o *X_X_TURN* - possible turn values for two successive residues
     where neither is a proline


*/

char *GLY_COIL_MS[93]; 
char *COIL_MS[94];
char *PRO_COIL_MS[18];
char *TURN_MS[54][2];
char *PRO_PRO_TURN[4][2]; 
char *PRO_GLY_TURN[3][2]; 
char *PRO_X_TURN[15][2]; 
char *GLY_PRO_TURN[2][2];
char *X_PRO_TURN[5][2];
char *PII_MS[2];
char *HELIX_MS[6];
char *STRAND_MS[9];
char *PRO_STRAND_MS[1];

