#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double linus_sy(double psi)
{
    double sy = psi + 180.0;
    if(sy > 180.0)
        sy = sy - 360.0;
    else if (sy < -180.0)
        sy = sy + 360.0;

    return sy;
}

