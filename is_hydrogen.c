#include <stdio.h>
#include <string.h>
int is_hydrogen(char * name)
{
    int is_H =0;
    int i;
    char *hydrogen[] = {" H  "," HA "," HB "," HD1"," HD2"," HE "," HE1"," HE2"," HE3"," HG "," HG1"," HH "," HH2"," HZ "," HZ2"," HZ3","1H  ","1HA ","1HB ","1HD ","1HD1","1HD2","1HE ","1HE2","1HG ","1HG1","1HG2","1HH1","1HH2","1HH3","1HW1","1HW2","1HZ ","2H  ","2HA ","2HB ","2HD ","2HD1","2HD2","2HE ","2HE2","2HG ","2HG1","2HG2","2HH1","2HH2","2HH3","2HW1","2HW2","2HZ ","3H  ","3HB ","3HD1","3HD2","3HE ","3HG1","3HG2","3HH3","3HZ "}; 
    //printf("-%s-\n",name);
    
    for(i=0;i<59;i++)
    {
       
        if(!strcmp(hydrogen[i],name))
        {
            //printf("Found match!\n");
            is_H=1;
            break;
        }
        
    }    
    return is_H;
}
