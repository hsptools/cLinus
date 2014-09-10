#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Since for loops have ctr < end, increased end index by 1. For eg end for DUM is 4 while ACE is 7

int get_res_range(char *res_name, int * start, int * end)
{
   if (!strcmp(res_name, "DUM"))
   {
     *start = 0;
      *end = 5;
   }
   else if (!strcmp(res_name, "ACE"))
   {
      *start = 5;
      *end = 8;
   }
   else if (!strcmp(res_name, "NH2"))
   {
      *start = 8;
      *end = 9;
   }
   else if (!strcmp(res_name, "NME"))
   {
      *start = 9;
      *end = 12;
   }
   else if (!strcmp(res_name, "ALA"))
   {
      *start = 12;
      *end = 19;
   }
   else if (!strcmp(res_name, "ARG"))
   {
      *start = 19;
      *end = 32;
   }
   else if (!strcmp(res_name, "ASN"))
   {
      *start = 32;
      *end = 42;
   }
   else if (!strcmp(res_name, "ASP"))
   {
      *start = 42;
      *end = 52;
   }
   else if (!strcmp(res_name, "CYS"))
   {
      *start = 52;
      *end = 60;
   }
   else if (!strcmp(res_name, "GLN"))
   {
      *start = 60;
      *end = 71;
   }
   else if (!strcmp(res_name, "GLU"))
   {
      *start = 71;
      *end = 82;
   }
   else if (!strcmp(res_name, "GLY"))
   {
      *start = 82;
      *end = 89;
   }
   else if (!strcmp(res_name, "HIS"))
   {
      *start = 89;
      *end = 101;
   }
   else if (!strcmp(res_name, "ILE"))
   {
      *start = 101;
      *end = 111;
   }
   else if (!strcmp(res_name, "LEU"))
   {
      *start = 111;
      *end = 121;
   }
   else if (!strcmp(res_name, "LYS"))
   {
      *start = 121;
      *end = 132;
   }
   else if (!strcmp(res_name, "MET"))
   {
      *start = 132;
      *end = 142;
   }
   else if (!strcmp(res_name, "PHE"))
   {
      *start = 142;
      *end = 155;
   }
   else if (!strcmp(res_name, "PRO"))
   {
      *start = 155;
      *end = 163;
   }
   else if (!strcmp(res_name, "SER"))
   {
      *start = 163;
      *end = 171;
   }
   else if (!strcmp(res_name, "THR"))
   {
      *start = 171;
      *end = 180;
   }
   else if (!strcmp(res_name, "TRP"))
   {
      *start = 180;
      *end = 196;
   }
   else if (!strcmp(res_name, "TYR"))
   {
      *start = 196;
      *end = 210;
   }
   else if (!strcmp(res_name, "VAL"))
   {
      *start = 210;
      *end = 219;
   }
   else
   {
      printf("Residue name not found!\n");
      exit(1);
   }
}


