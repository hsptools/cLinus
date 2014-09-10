#include <stdio.h>
#include <string.h>
#include "linus.h"

/*
 Methods of this class

        o *get_logfile* - returns a reference to the LOG file

        o *get_out_file* - returns a reference to the file used for saving
        structures during a simulation

        o *get_min_file* - returns a reference to the file used for saving
        the best energy structure found during a simulation

        o *get_rc_file* - returns a reference to the file used for saving
        the mesostate conformation strings found during a simulation

        o *get_swt_file* - returns *name* of the file used for saving
        the observed distribution of the residues in the various
        secondary structure states

        o *get_ps_file* - returns *name* of the file to which the
        postscript plot of the the observed distribution of the residues in
        the various secondary structure states is written

        o *create_file_name* - given an extension, creates the *name* of
        the file

*/


void create_linusfile(LinusFile *lf,char *pdbfilename)
{
    strcpy(lf->input,pdbfilename);
    
    char *ext = strtok(pdbfilename,".pdb");
    strcpy(lf->basename,ext);
    
    strcpy(lf->log_filename,lf->basename);
    strcat(lf->log_filename,".LOG");
    
    strcpy(lf->out_filename,lf->basename);
    strcat(lf->out_filename,".out");
    
    strcpy(lf->min_filename,lf->basename);
    strcat(lf->min_filename,".min");
    
    strcpy(lf->rc_filename,lf->basename);
    strcat(lf->rc_filename,".rc");
    
    strcpy(lf->swt_filename,lf->basename);
    strcat(lf->swt_filename,".swt");
    
    strcpy(lf->msd_filename,lf->basename);
    strcat(lf->msd_filename,".msd");
    
    strcpy(lf->pms_filename,lf->basename);
    strcat(lf->pms_filename,".pms");
    
    strcpy(lf->ps_filename,lf->basename);
    strcat(lf->ps_filename,".ps");
}



