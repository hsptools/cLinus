#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linus.h"

void print_conf(char *fmt, LinusProtein *prot, FILE *ofile,int num_steric_hbs)
{
    double phi1=prot->atoms[prot->phi_atoms[1]].tp_torsion;
    double phi2=prot->atoms[prot->phi_atoms[2]].tp_torsion;
    double phi3=prot->atoms[prot->phi_atoms[3]].tp_torsion;
    double phi4=prot->atoms[prot->phi_atoms[4]].tp_torsion;
       
    double psi1=out_sy(prot->atoms[prot->psi_atoms[1]].tp_torsion);
    double psi2=out_sy(prot->atoms[prot->psi_atoms[2]].tp_torsion);
    double psi3=out_sy(prot->atoms[prot->psi_atoms[3]].tp_torsion);
    double psi4=out_sy(prot->atoms[prot->psi_atoms[4]].tp_torsion);
    
    double ome1=prot->atoms[prot->omega_atoms[1]].tp_torsion;
    double ome2=prot->atoms[prot->omega_atoms[2]].tp_torsion;
    double ome3=prot->atoms[prot->omega_atoms[3]].tp_torsion;
    double ome4=prot->atoms[prot->omega_atoms[4]].tp_torsion;
    
    
    //printf(fmt,phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3,  phi4, psi4, ome4, num_steric_hbs);    
    fprintf(ofile,fmt,phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3,  phi4, psi4, ome4, num_steric_hbs);
}

