#include "pileup_stat_lib.h"
    
int main(int argc, char** argv)
{
    //cmdline args
    //parameters for mutation calling
    if(argc<2){
        printf("ERROR please provide  baseq limit,\
               and a list of sample names \n");
        exit(1);
    }
    
    //int cov_limit=(int) strtol(argv[1],NULL,10);
    int baseq_limit=(int) strtol(argv[1],NULL,10);

    //sample names
    int n_sample_names=argc-2;
    char** sample_names= (char**) malloc(n_sample_names * sizeof(char*));
    int i=0;
    for(i=0;i<n_sample_names;i++) sample_names[i]=argv[2+i];
    
    //variables for reading a line
    char* line = NULL;
    size_t len = 0;
    ssize_t line_size;
    
    //print header
    for(i=0;i<n_sample_names-1;i++){
        printf("%s_A\t%s_C\t%s_G\t%s_T\t",sample_names[i],
               sample_names[i],sample_names[i],sample_names[i]);
    }
    printf("%s_A\t%s_C\t%s_G\t%s_T\n",sample_names[i],
           sample_names[i],sample_names[i],sample_names[i]);
    
    //loop over input lines
    //FILE* test_f = fopen("/Users/ribli/pileup_stat_lib/dbg.pup","r");
    //while ((line_size = getline(&line, &len, test_f)) != -1) {
    while ((line_size = getline(&line, &len, stdin)) != -1) {
        //the pileup line structure for the line being read
        struct mplp my_mplp;
        init_mplp(&my_mplp);
        
        //build the struct from input line
        process_mplp_input_line(&my_mplp,line,line_size,baseq_limit,sample_names,n_sample_names);
        
        //print the freq
        print_counts(&my_mplp);
        
        //free the memory allocated by the struct
        free_mplp(&my_mplp);
    }
    //free resources
    free(line);
    return 0;
}
