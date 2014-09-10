/*
 """Class to create files that are output by LINUS.  All filenames
    use the base name of the input file and then add the appropriate
    extensions to create new files.

   
    """
*/

typedef struct
{
    char input[200];
    char basename[200];
    char log_filename[200];
    //FILE *log_file;
    char out_filename[200];
    //FILE *out_file;
    char min_filename[200];
    //FILE *min_file;
    char rc_filename[200];    
    //FILE *rc_file;
    char swt_filename[200];    
    //FILE *swt_file;
    char msd_filename[200];    
    //FILE *msd_file;
    char pms_filename[200];    
    //FILE *pms_file;
    char ps_filename[200];    
    //FILE *ps_file;
}LinusFile;

