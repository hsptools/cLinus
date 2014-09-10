#include <stdio.h>
//#include "Atom3d.h"
#include "linus.h"
void restore_coords(Atom3d *atom)
{
   //printf("OLD %lf, %lf, %lf\n",atom->fp_distance, atom->sp_angle, atom->tp_torsion);
   atom->fp_distance = atom->fp_distance_tmp;
   atom->sp_angle = atom->sp_angle_tmp;
   atom->tp_torsion = atom->tp_torsion_tmp;
   //printf("NEW %lf, %lf, %lf\n",atom->fp_distance, atom->sp_angle, atom->tp_torsion);

}

void restore_coords_range(Atom3d *atoms, int start, int end) // start=None, end=None)
{
   // Save internal coordinates following a conformation change
   int i;
   for (i = start; i < end; i++)
   {
       restore_coords(&atoms[i]);

   }
}

