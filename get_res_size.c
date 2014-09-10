#include <stdio.h>
#include <stdlib.h>

int get_res_size(char *res_name)
{
   if (!strcmp(res_name, "DUM"))
   {
      //numres[0] = 0;
      //numres[1] = 4;
      return 5;
   }
   else if (!strcmp(res_name, "ACE"))
   {
      //numres[0] = 5;
      //numres[1] = 7;
      return 3;
   }
   else if (!strcmp(res_name, "NH2"))
   {
      //numres[0] = 8;
      //numres[1] = 8;
      return 1;
   }
   else if (!strcmp(res_name, "NME"))
   {
      //numres[0] = 9;
      //numres[1] = 11;
      return 3;
   }
   else if (!strcmp(res_name, "ALA"))
   {
      //numres[0] = 12;
      //numres[1] = 18;
      return 7;
   }
   else if (!strcmp(res_name, "ARG"))
   {
      //numres[0] = 19;
      //numres[1] = 31;
      return 13;
   }
   else if (!strcmp(res_name, "ASN"))
   {
      //numres[0] = 32;
      //numres[1] = 41;
      return 10;
   }
   else if (!strcmp(res_name, "ASP"))
   {
      //numres[0] = 42;
      //numres[1] = 51;
      return 10;
   }
   else if (!strcmp(res_name, "CYS"))
   {
      //numres[0] = 52;
      //numres[1] = 59;
      return 8;
   }
   else if (!strcmp(res_name, "GLN"))
   {
      //numres[0] = 60;
      //numres[1] = 70;
      return 11;
   }
   else if (!strcmp(res_name, "GLU"))
   {
      //numres[0] = 71;
      //numres[1] = 81;
      return 11;
   }
   else if (!strcmp(res_name, "GLY"))
   {
      //numres[0] = 82;
      //numres[1] = 88;
      return 7;
   }
   else if (!strcmp(res_name, "HIS"))
   {
      //numres[0] = 89;
      //numres[1] = 100;
      return 12;
   }
   else if (!strcmp(res_name, "ILE"))
   {
      //numres[0] = 101;
      //numres[1] = 110;
      return 10;
   }
   else if (!strcmp(res_name, "LEU"))
   {
      //numres[0] = 111;
      //numres[1] = 120;
      return 10;
   }
   else if (!strcmp(res_name, "LYS"))
   {
      //numres[0] = 121;
      //numres[1] = 131;
      return 11;
   }
   else if (!strcmp(res_name, "MET"))
   {
      //numres[0] = 132;
      //numres[1] = 141;
      return 10;
   }
   else if (!strcmp(res_name, "PHE"))
   {
      //numres[0] = 142;
      //numres[1] = 154;
      return 13;
   }
   else if (!strcmp(res_name, "PRO"))
   {
      //numres[0] = 155;
      //numres[1] = 162;
      return 8;
   }
   else if (!strcmp(res_name, "SER"))
   {
      //numres[0] = 163;
      //numres[1] = 170;
      return 8;
   }
   else if (!strcmp(res_name, "THR"))
   {
      //numres[0] = 171;
      //numres[1] = 179;
      return 9;
   }
   else if (!strcmp(res_name, "TRP"))
   {
      //numres[0] = 180;
      //numres[1] = 195;
      return 16;
   }
   else if (!strcmp(res_name, "TYR"))
   {
      //numres[0] = 196;
      //numres[1] = 209;
      return 14;
   }
   else if (!strcmp(res_name, "VAL"))
   {
      //numres[0] = 210;
      //numres[1] = 218;
      return 9;
   }
   else
   {
      printf("Residue name not found!\n");
      exit(1);
   }
}
