#include "pileup_stat_lib.h"
#include <math.h>
#include <math.h>

    
#define NPARAMS 1
    
int main(int argc, char** argv)
{
    //cmdline args
    if(argc!= NPARAMS + 1){
        printf("ERROR please provide %d args \n", NPARAMS); 
        exit(1);
    }
    
    //parameters 
    
    //varaiables for reading a line
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    
    //the pileup structure
    struct Mpileup_line my_pup_line;
    init_mpileup_line(&my_pup_line);
    
    //loop over input lines
    while ((read = getline(&line, &len, stdin)) != -1) {
        //read mpileup line
        get_mpileup_line(&my_pup_line,line,read);
       
        //count the bases
        count_bases_all_sample(&my_pup_line,baseq_limit);
        //calcluate freqs
        calculate_base_freqs_all_sample(&my_pup_line);
        
    }

    //free resources
    free_mpileup_line(&my_pup_line);
    free(line);
    return 0;
}
