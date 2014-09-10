#include "linusSetConf.h"
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



char *GLY_COIL_MS[93] = {"Aa","Ab","Ak","Al",
               "Ba","Bb","Bf","Bg","Bh","Bj","Bk","Bl",
               "Ca","Cb","Cc","Cd","Ce","Cf","Cg","Ch","Cj","Ck","Cl",
               "Da","Db","Dc","Dd","De","Df","Dg","Dh","Di","Dj","Dk","Dl",
               "Ea","Eb","Ed","Ee","Ef","Eg","Ej","Ek","El",
               "Fd","Fe","Ff","Fk",
               "Hb","Hc","Hg","Hh","Hi",
               "Ia","Ib","Ic","Ie","If","Ig","Ih","Ii","Ik","Il",
               "Ja","Jb","Jc","Jd","Je","Jf","Jg","Jh","Ji","Jj","Jk","Jl",
               "Ka","Kb","Kc","Ke","Kf","Kg","Kh","Ki","Kj","Kk","Kl",
               "La","Lb","Le","Lf","Lg","Lk","Ll"};

char *COIL_MS[94] = {"Aa","Aj","Ak","Al",
           "Ba","Bb","Bc","Bd","Be","Bf","Bg","Bh","Bi","Bj","Bk","Bl",
           "Ca","Cb","Cc","Cd","Ce","Cf","Cg","Ch","Ci","Cj","Ck","Cl",
           "Da","Db","Dc","Dd","De","Df","Dg","Dh","Di","Dj","Dk","Dl",
           "Ea","Ed","Ee","Ef","Eg","Ej","Ek","El",
           "Fd","Fe","Ff","Fj","Fk","Fl",
           "Hb","Hc","Hg","Hh","Hi",
           "Ia","Ib","Ic","Id","Ie","If","Ig","Ih","Ii","Il",
           "Jd","Je","Jf","Jg","Jh","Ji","Jl",
           "Dl","Dl","Dl","Dl","Dl","Dl","Dl","Dl","Dl",
           "Di","Di","Di","Di","Di","Di","Di","Di","Di"};

char *PRO_COIL_MS[18] = {"De","Df","Dg","Dh","Di","Dj","Dk","Dl",
               "Ed","Ee","Ef","Ef","Eg","Eh","Ei","Ej","Ek","El"};
               
char *PII_MS[2] = {"Ek","El"};

char *HELIX_MS[6] = {"Ee","Ee","Ee","De","Df","Ef"};

char *STRAND_MS[9] = {"Ck","Ck","Ck","Cj","Cl","Bk","Bl","Dk","Dj"};

char *PRO_STRAND_MS[1] = {"Ek"}; 

// does not form many turns - a few type 3
// make most choices PII

char *PRO_PRO_TURN [4][2] = {
                {"Ee","Ee"},
                {"Ek","Ek"},
                {"Ek","El"},
                {"El","Ek"},
               };
char *PRO_GLY_TURN [3][2]= {
                {"Ek","Jf"},
                {"Ek","Ig"},
                {"Ek","Jg"},
               };

char  *PRO_X_TURN[15][2] = {
              {"Ef","Ef"},{"Ef","Ef"},{"Ef","Ef"},
              {"Ef","Dg"},{"Ef","Dg"},
              {"Ee","Dg"},{"Ee","Dg"},
              {"Ee","Ee"},{"Ee","Ee"},
              {"Ef","Cg"},{"Ef","Cg"},
              {"Ee","Df"},{"Ee","Df"},
              {"Ek","Jf"},
              {"Ek","Ig"}
             };

char  *GLY_PRO_TURN [2][2]= {
                {"Ib","Ee"},
                {"Ib","Ef"}
               };

// does not form many turns - upstream of PII is left udder
char  *X_PRO_TURN[5][2] = {
              {"Ee","Ef"},
              {"Bi","Ek"},
              {"Ci","Ek"},
              {"Dl","Ek"},
              {"El","Ek"}
             };




char *TURN_MS[54][2] = {{"Ef","Df"},{"Ef","Df"},{"Ef","Df"},
           {"Ee","Ef"},{"Ee","Ef"},{"Ee","Ef"},
           {"Ef","Ef"},{"Ef","Ef"},{"Ef","Ef"},
           {"Ef","Dg"},{"Ef","Dg"},
           {"Ee","Dg"},{"Ee","Dg"},
           {"Ee","Ee"},{"Ee","Ee"},
           {"Ef","Cg"},{"Ef","Cg"},
           {"Ee","Df"},{"Ee","Df"},
           {"Ek","Jf"},
           {"Ek","Ig"},
           {"Ef","Ee"},
           {"Ek","Jg"},
           {"Ee","Cg"},
           {"Df","Df"},
           {"Ef","Cf"},
           {"Dg","Df"},
           {"Df","Dg"},
           {"Ih","Ig"},
           {"Ef","De"},
           {"Ek","Ih"},
           {"Dg","Cg"},
           {"Df","Cg"},
           {"Ib","Dg"},
           {"Df","Ee"},
           {"Fe","Ef"},
           {"Ib","Ef"},
           {"Df","Ef"},
           {"Ih","Jf"},
           {"Ih","Jg"},
           {"Ig","Ig"},
           {"Ef","Ch"},
           {"Dg","Ee"},
           {"Dg","Ef"},
           {"Ee","Eg"},
           {"Ih","Ih"},
           {"Ee","De"},
           {"Ig","Jg"},
           {"Ig","Jf"},
           {"Ib","Df"},
           {"Ic","Df"},
           {"Ic","Dg"},
           {"Dg","Dg"},
           {"Ig","Ih"},
          };
