#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linusSecStr.h"
#include "linus.h"

//res_rc returns mesostate code given phi/psi

void res_rc(double r1, double r2, double r3,char *code)
{
    //printf("%lf %lf %lf\n",r1,r2,r3);
    if(r3 ==0.0) r3 = 180.0;
    if(fabs(r3) <= 90.0)
    {
        strcpy(code,"Xx");
    }
    else if(r1 > 180.0 || r2 > 180.0)
    {
        strcpy(code,"??");
    }
    printf("r1 %lf %d\n",r1/30.0,round1(r1/30.0));
    int ir1 = (int)(round1(r1/30.0)) *30;
        
    //special case of psi equals exactly 180
    if((int)r2 == 180)
    {
        r2 = 179.0;
    }
    else if((int)r2 == -180)
    {
        r2=-179.0;
    }
    printf("r2 %lf %d\n",((r2+15.0)/30.0),round1(((r2+15.0)/30.0)));
    int ir2 = -15 + (int)(round1((r2 + 15.0)  /30.0 )) * 30;
    
    if(ir1 == 80)
    {
        ir1 = -180;
    }
    printf("ir1 %d ir2 %d\n",ir1,ir2);
    get_rc_code(ir1,ir2,code);
    //printf("Code is %s\n",code);

}

int round1(double x)
{
    if(x > 1)
    {
        double x1=x*10.0;
        double y = (int)x1%10;
        if(y>=5) return (int)x+1.0;
        else return (int) x;
    }
    else if (x>0 && x<1)
    {
    
    
    
    }
    else
    {
        x = fabs(x);
        double x1=x*10.0;
        double y = (int)x1%10;
        if(y>=5) return (int)-(x+1.0);
        else return (int)-x;
    
    }
}
