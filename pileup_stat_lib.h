#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <ctype.h>
    
///////////////////////////////////////////////////////////////////////////
// Constants 
////////////////////////////////////////////////////////////////////////////    
   
//Size limitations
    
//maximum sample limit, change it if you have more samples
#define MAXSAMPLE 1024
//maximum length of an indel
#define MAX_INDEL_LEN 128 
   
    
//values for easier usage
    
//counts are stored in an array
//element indices are defined here
#define COV_IDX 0
#define REF_IDX 1 
#define A_IDX 2 
#define C_IDX 3 
#define G_IDX 4 
#define T_IDX 5 
#define DEL_IDX 6 
#define INS_START_IDX 7 
#define DEL_START_IDX  8 
#define READ_START_IDX 9 
#define READ_END_IDX 10 
#define MAX_IDX 11
   

// 0 coverage samples also have some kind of "frequency"
// should always check for this when working with the frequencies
#define ZERO_COV_FREQ -9999

    
///////////////////////////////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////////////////////////////////    
    
/*
    the mpileup struct
*/
struct mplp
{
    //sample names for easier output
    char** sample_names;
    int n_sample_names;
    
    //raw data
    char* raw_line;
    
    //position data
    char* chrom;
    int pos;
    char ref_nuq;
    
    //sample data
    int n_samples;
    
    //sample level raw data
    int raw_cov[MAXSAMPLE];
    char* raw_bases[MAXSAMPLE];
    char* raw_quals[MAXSAMPLE];
    
    //more structured data
    int counts[MAXSAMPLE][MAX_IDX];
    double freqs[MAXSAMPLE][MAX_IDX];
    char** ins_bases[MAXSAMPLE];
    char** del_bases[MAXSAMPLE];
    
    //mutation data 
    char mut_type[4];
    char mut_base;
    double mut_score;
    char mut_indel[MAX_INDEL_LEN];
    int mut_sample_idx;
    double mut_freq;
    double cleanliness;
    
};


////////////////////////////////////////////////////////////////////////////
// initializing, deleting, and copying pileup struct
////////////////////////////////////////////////////////////////////////////

/*
    initializes pointers with NULL
*/
int init_mplp(struct mplp* my_mplp);

/*
    frees all malloced objects in mplp struct, 
    and initializes them to NULL
*/
int free_mplp(struct mplp* my_mplp);

/*
    deep copy mpileup line
*/
int copy_mplp(struct mplp* target ,struct mplp* from);


////////////////////////////////////////////////////////////////////////////
// formatted printing of pileup struct
////////////////////////////////////////////////////////////////////////////

/*
    printf formatted representation of mpileup line
*/
int print_mplp(struct mplp* my_mplp);


////////////////////////////////////////////////////////////////////////////
// process input mpileup line from samtools mpileup command output
////////////////////////////////////////////////////////////////////////////

/*
      process input mpileup line from samtools mpileup command output
*/
int process_mplp_input_line(struct mplp* my_mplp,char* line, ssize_t line_size,
                            int baseq_limit,char** sample_names,int n_sample_names);

////////////////////////////////////////////////////////////////////////////
// read pileup struct from mpileup line
////////////////////////////////////////////////////////////////////////////


/*
    gets next entry from tab separated input string as null terminated char*
*/
int get_next_entry(char* line, ssize_t line_size, ssize_t* pointer, char** result);

/*
    gets mpileup line from the char* line
*/
int get_mplp(struct mplp* my_mplp,char* line, ssize_t read);


////////////////////////////////////////////////////////////////////////////
// Count bases in mplieup struct
////////////////////////////////////////////////////////////////////////////

/*
    counts bases in all samples
*/
int count_bases_all_samples(struct mplp* my_mplp,int baseq_lim);


/*
    counts bases in one sample
*/
int count_bases(char* bases,char* quals,int* counts,char ref_base,int baseq_lim);

/* 
    parse a base from the bases and quals
*/
int handle_base(char* bases,char* quals,int* counts, int* filtered_cov,int* base_ptr,int* qual_ptr,int baseq_lim);


/* 
    parse a deletion from the bases and quals
*/
int handle_deletion(char* bases,int* del_count,int* base_ptr, char qual,int baseq_lim);

/* 
    parse an insertion from the bases and quals
*/
int handle_insertion(char* bases,int* ins_count,int* base_ptr, char qual,int baseq_lim);


////////////////////////////////////////////////////////////////////////////
// Calculate base + ins + del frequencies
////////////////////////////////////////////////////////////////////////////


/*
    calculate base freqs in all samples
*/
int calculate_freqs_all_samples(struct mplp* my_mplp);


/*
   calculate freqs in a sample
*/
int calculate_freqs(double* freqs,int* counts);


////////////////////////////////////////////////////////////////////////////
// Collect indels
////////////////////////////////////////////////////////////////////////////

/*
    collect indels in all samples
*/
int collect_indels_all_samples(struct mplp* my_mplp,int baseq_lim);

/*
    collect the inserted, and deleted bases
*/
int collect_indels(char* bases,char* quals, char*** ins_bases, int ins_count, 
                   char*** del_bases,int del_count, int baseq_lim);

/*
    free memory of indel bases
*/
int free_indel_bases(char*** ins_bases,char*** del_bases);


////////////////////////////////////////////////////////////////////////////
// Print stuff 
////////////////////////////////////////////////////////////////////////////
