//HbList
//#include "Acps.h"
//Acps
typedef struct
{
    int acp;
    int acp1;
    double hbdist;
    double hbdmax;
    double hbtor;
    double hbene;
}Acps;


typedef struct
{
    int donor;
    int da1;
    int da2;
    Acps *acps;
    int acps_size;
}HbList;

