#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linus.h"


int get_hexname(double dh1, double dh2, double dh3)
{
    int step = 60;
    int hex_num = 0;
    int hexname;
    int a,b,c;
    for (a=-180;a<180;a+=step)
    {
        for(b=-180;b<180;b+=step)
        {
            for(c=-180;c<180;c+=step)
            {
                hex_num+=1;
                if(dh1 >= a && dh1 <= a+step)
                {
                    if(dh2 >= b && dh2 <= b+step)
                    {
                        if(dh3 >= c && dh3 <= c+step)
                        {
                            hexname = hex_num;
                        }   
                    }   
                }
            }
        }
    }
    return hexname;
}
