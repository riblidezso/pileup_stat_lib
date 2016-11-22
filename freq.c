#include "pileup_stat_lib.h"
    
int main(int argc, char** argv)
{
    //cmdline args
    //parameters for mutation calling
    if(argc<1){
        printf("ERROR please provide  list of sample names \n");
        exit(1);
    }
    
    int cov_limit=(int) strtol(argv[3],NULL,10);
    int baseq_limit=(int) strtol(argv[4],NULL,10);

    //sample names
    int n_sample_names=argc-1;
    char** sample_names= (char**) malloc(n_sample_names * sizeof(char*));
    int i=0;
    for(i=0;i<n_sample_names;i++) sample_names[i]=argv[1+i];
    
    //variables for reading a line
    char* line = NULL;
    size_t len = 0;
    ssize_t line_size;
    
    //print header
    //printf("#sample_name\tchr\tpos\ttype\tscore\tref\tmut\tcov\tmut_freq\tcleanliness\n");
    
    //loop over input lines
    //FILE* test_f = fopen("/Users/ribli/unique_mutation/data2/dbg.input","r");
    //while ((line_size = getline(&line, &len, test_f)) != -1) {
    while ((line_size = getline(&line, &len, stdin)) != -1) {
        //the pileup line structure for the line being read
        struct mplp my_mplp;
        init_mplp(&my_mplp);
        
        //build the struct from input line
        process_mplp_input_line(&my_mplp,line,line_size,baseq_limit,sample_names,n_sample_names);
        
        //free the memory allocated by the struct
        free_mplp(&my_mplp);
    }
    //free resources
    free(line);
    return 0;
}