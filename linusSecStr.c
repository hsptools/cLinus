#include <stdio.h>
#include "linusSecStr.h"

int rc_phi[144] = {-180,-180,-180,-180,-180,-180,-180,-180,-180,-180,-180,-180,-150,-150,-150,-150,-150,-150,-150,-150,-150,-150,-150,-150,-120,-120,-120,-120,-120,-120,-120,-120,-120,-120,-120,-120,-90,-90,-90,-90,-90,-90,-90,-90,-90,-90,-90,-90,-60,-60,-60,-60,-60,-60,-60,-60,-60,-60,-60,-60,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,0,0,0,0,0,0,0,0,0,0,0,0,30,30,30,30,30,30,30,30,30,30,30,30,60,60,60,60,60,60,60,60,60,60,60,60,90,90,90,90,90,90,90,90,90,90,90,90,120,120,120,120,120,120,120,120,120,120,120,120,150,150,150,150,150,150,150,150,150,150,150,150};
int rc_psi[144] = {-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165,-165,-135,-105,-75,-45,-15,15,45,75,105,135,165};
char *rc_code[144]={"Aa", "Ab", "Ac", "Ad", "Ae", "Af", "Ag", "Ah", "Ai", "Aj", "Ak", "Al", "Ba", "Bb", "Bc", "Bd", "Be", "Bf", "Bg", "Bh", "Bi", "Bj", "Bk", "Bl", "Ca", "Cb", "Cc", "Cd", "Ce", "Cf", "Cg", "Ch", "Ci", "Cj", "Ck", "Cl", "Da", "Db", "Dc", "Dd", "De", "Df", "Dg", "Dh", "Di", "Dj", "Dk", "Dl", "Ea", "Eb", "Ec", "Ed", "Ee", "Ef", "Eg", "Eh", "Ei", "Ej", "Ek", "El", "Fa", "Fb", "Fc", "Fd", "Fe", "Ff", "Fg", "Fh", "Fi", "Fj", "Fk", "Fl", "Ga", "Gb", "Gc", "Gd", "Ge", "Gf", "Gg", "Gh", "Gi", "Gj", "Gk", "Gl", "Ha", "Hb", "Hc", "Hd", "He", "Hf", "Hg", "Hh", "Hi", "Hj", "Hk", "Hl", "Ia", "Ib", "Ic", "Id", "Ie", "If", "Ig", "Ih", "Ii", "Ij", "Ik", "Il", "Ja", "Jb", "Jc", "Jd", "Je", "Jf", "Jg", "Jh", "Ji", "Jj", "Jk", "Jl", "Ka", "Kb", "Kc", "Kd", "Ke", "Kf", "Kg", "Kh", "Ki", "Kj", "Kk", "Kl", "La", "Lb", "Lc", "Ld", "Le", "Lf", "Lg", "Lh", "Li", "Lj", "Lk", "Ll"};

char *PII_meso[4] = {"Dk","Dl","Ek", "El"};

int is_PII(char * meso)
{
    int i;
    for(i=0;i<4;i++)
    {
        if(!strcmp(PII_meso[i],meso))
        {
            return 1;
        }
    }
    return 0;
}

char *h_meso[6] = {"De","Df","Ed","Ee","Ef","Fe"};

int is_helix(char * meso)
{
    int i;
    for(i=0;i<6;i++)
    {
        if(!strcmp(h_meso[i],meso))
        {
            return 1;
        }
    }
    return 0;
}

char *s_meso[8]={"Bj","Bk","Bl","Cj","Ck","Cl","Dj","Dk"};

int is_strand(char * meso)
{
    int i;
    for(i=0;i<8;i++)
    {
        if(!strcmp(s_meso[i],meso))
        {
            return 1;
        }
    }
    return 0;
}

char *t_meso[144] = {
         "EfDf",  //5226  relative number of occurance in PDB
         "EeEf",  //4593 
         "EfEf",  //4061 
         "EfDg",  //3883 
         "EeDg",  //2118 
         "EeEe",  //1950 
         "EfCg",  //1932 
         "EeDf",  //1785 
         "EkJf",  //1577 
         "EkIg",  //1106 
         "EfEe",  // 995 
         "EkJg",  // 760 
         "EeCg",  // 553 
         "DfDf",  // 479 
         "EfCf",  // 395 
         "DgDf",  // 332 
         "DfDg",  // 330 
         "IhIg",  // 310 
         "EfDe",  // 309 
         "EkIh",  // 298 
         "DgCg",  // 275 
         "DfCg",  // 267 
         "IbDg",  // 266 
         "DfEe",  // 260 
         "FeEf",  // 250 
         "IbEf",  // 249 
         "DfEf",  // 219 
         "IhJf",  // 216 
         "IhJg",  // 213 
         "IgIg",  // 207 
         "EfCh",  // 188 
         "DgEe",  // 180 
         "DgEf",  // 176 
         "EeEg",  // 172 
         "IhIh",  // 153 
         "EeDe",  // 150 
         "IgJg",  // 147 
         "EkKf",  // 147 
         "EeCh",  // 147 
         "IbDf",  // 131 
         "DgDg",  // 128 
         "EgDf",  // 127 
         "FeDg",  // 114 
         "ElIg",  // 111 
         "IgIh",  // 107 
         "DfDe",  // 107 
         "EjIg",  // 101 
         "EeCf",  // 100 
         "DfCh",  //  94 
         "DgCf",  //  91 
         "DfCf",  //  91 
         "DeEe",  //  91 
         "DkIh",  //  88 
         "FeDf",  //  79 
         "EkIf",  //  78 
         "EeDh",  //  76 
         "DgCh",  //  74 
         "IgJf",  //  71 
         "EjJg",  //  71 
         "FeEe",  //  69 
         "DlIh",  //  66 
         "EgCg",  //  65 
         "ElIh",  //  62 
         "EjJf",  //  62 
         "FeCg",  //  59 
         "DlIg",  //  56 
         "IbCg",  //  54 
         "EfEg",  //  54 
         "EkJe",  //  53 
         "FkJf",  //  52 
         "ElJg",  //  51 
         "DgDe",  //  49 
         "DlJg",  //  46 
         "EgCf",  //  45 
         "IaEf",  //  40 
         "FkIg",  //  39 
         "JaEf",  //  38 
         "EjIh",  //  38 
         "EgEf",  //  38 
         "DkJg",  //  36 
         "DeEf",  //  34 
         "EeCi",  //  31 
         "JgIh",  //  29 
         "IcEf",  //  29 
         "EkKe",  //  29 
         "DkIg",  //  29 
         "IbEe",  //  27 
         "EgDg",  //  27 
         "EeFe",  //  27 
         "EjKf",  //  26 
         "IaDf",  //  25 
         "HhIg",  //  24 
         "HbDg",  //  24 
         "ElJf",  //  24 
         "EfDh",  //  24 
         "IcDf",  //  23 
         "EfBh",  //  23 
         "IcDg",  //  22 
         "IcCg",  //  22 
         "FkJg",  //  21 
         "FeCh",  //  21 
         "IgKf",  //  20 
         "FdDg",  //  20 
         "EkHh",  //  20 
         "DfDh",  //  20 
         "DgBh",  //  19 
         "DfBh",  //  19 
         "DeDf",  //  19 
         "DfFe",  //  18 
         "EfFe",  //  17 
         "EgEe",  //  16 
         "EgDe",  //  16 
         "DkJf",  //  16 
         "JgJg",  //  15 
         "IbEg",  //  15 
         "IbCh",  //  15 
         "EfBg",  //  15 
         "DgCe",  //  15 
         "JlEf",  //  14 
         "CgCg",  //  14 
         "HhJf",  //  13 
         "EeBi",  //  13 
         "DfBi",  //  13 
         "IhIf",  //  12 
         "FeEg",  //  12 
         "FdEf",  //  12 
         "EdEf",  //  12 
         "DlJf",  //  12 
         "DhCg",  //  12 
         "JgIg",  //  11 
         "IeBg",  //  11 
         "FjIg",  //  11 
         "FdCh",  //  11 
         "EdEe",  //  11 
         "JfIh",  //  10 
         "JaEe",  //  10 
         "HhJg",  //  10 
         "HbEf",  //  10 
         "HbCh",  //  10 
         "FkIh",  //  10 
         "FjJf",  //  10 
         "ElJe",  //  10 
         "DhDf",  //  10 
         "CgDf"  //  10 
         };
         
         
int is_turn(char * meso)
{
    int i;
    for(i=0;i<144;i++)
    {
        //printf("  turn %s\n",t_meso[i]);
        if(!strcmp(t_meso[i],meso))
        {
            return 1;
        }
    }
    return 0;
}


