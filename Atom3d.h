typedef struct{
    double x;        /* atom X coordinate */
    double y;        /* atom Y coordinate */
    double z;        /* atom Z coordinate */
    double radius;        /* atom radius */
    char name[5];        /* atom name */
    int fp_ind; // p->atoms index
    int sp_ind;  //p->atoms index
    int tp_ind;  //p->atoms index
    char first_name[5];
    int fp; // //linusres offset
    char second_name[5];
    int sp; // //linusres offset
    char third_name[5];
    int tp; // //linusres offset
    double fp_distance;        /* distance of atom to it's first Parent */
    double sp_angle;        /* angle made by atom with first and second parents */
    double tp_torsion;        /* torsion angle between planes defined by atom, first parent, second parent AND first parent, second parent, third parent */
    double fp_distance_tmp;    /* distance of atom to it's first Parent */
    double sp_angle_tmp;    /* angle made by atom with first and second parents */
    double tp_torsion_tmp;    /* torsion angle between planes defined by atom, first parent, second parent AND first parent, second parent, third parent */
    int bsep0, bsep1;
    int resnum;
    double sasa;
    double sasarad;
    double sasapar;
}Atom3d;
