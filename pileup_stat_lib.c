#include "pileup_stat_lib.h"


////////////////////////////////////////////////////////////////////////////
// initialize
////////////////////////////////////////////////////////////////////////////
/*
    initializes all the pointers with NULL
        - always call this first after creating the struct,
            to avoid ERROR due to freeing/using memory garbage
*/
int init_mplp(struct mplp* my_mplp){
    int i;
    my_mplp->raw_line=NULL;
    my_mplp->chrom=NULL;
    my_mplp->pos=-42;
    my_mplp->n_samples=-42;
    my_mplp->ref_nuq='E';

    for(i=0;i<MAXSAMPLE;i++){
        my_mplp->raw_bases[i]=NULL;
        my_mplp->raw_quals[i]=NULL;
        my_mplp->ins_bases[i]=NULL;
        my_mplp->del_bases[i]=NULL;
    }


    strcpy(my_mplp->mut_type,"NOT\0");
    my_mplp->mut_base='E';
    my_mplp->mut_score=-42;
    for (i=0;i<MAX_INDEL_LEN;i++){my_mplp->mut_indel[i]='.';}
    my_mplp->mut_sample_idx=-42;
    my_mplp->mut_freq=-42;
    my_mplp->cleanliness=-42;

    return 0;
}

////////////////////////////////////////////////////////////////////////////
// free resources
////////////////////////////////////////////////////////////////////////////
/*
    frees all malloced objects in mplp struct
        - dont try to free uninitialized struct!
*/
int free_mplp(struct mplp* my_mplp){
    if( my_mplp->raw_line !=NULL) free(my_mplp->raw_line);
    if( my_mplp->chrom !=NULL) free(my_mplp->chrom);

    int i,j;
    for(i=0;i<my_mplp->n_samples;i++){
        if( my_mplp->raw_bases[i] !=NULL) free(my_mplp->raw_bases[i]);
        if( my_mplp->raw_quals[i] !=NULL) free(my_mplp->raw_quals[i]);

        if( my_mplp->ins_bases[i] !=NULL){
            for(j=0;j<my_mplp->counts[i][INS_START_IDX];j++){
                free(my_mplp->ins_bases[i][j]);
            }
            free(my_mplp->ins_bases[i]);
        }
        if( my_mplp->del_bases[i] !=NULL){
            for(j=0;j<my_mplp->counts[i][DEL_START_IDX];j++){
                free(my_mplp->del_bases[i][j]);
            }
            free(my_mplp->del_bases[i]);
        }
    }
    //initalize the pointers to null after freeing
    init_mplp(my_mplp);
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// deep copy pileup struct
////////////////////////////////////////////////////////////////////////////

/*
    deep copy mpileup line
*/
int copy_mplp(struct mplp* target ,struct mplp* from) {
    //free and initialize resources int the target
    free_mplp(target);

    target->sample_names=from->sample_names;
    target->n_sample_names=from->n_sample_names;

    //copy raw line and chromosome
    target->raw_line = (char*) malloc((strlen(from->raw_line)+1) * sizeof(char));
    strcpy( target->raw_line, from->raw_line);
    target->chrom = (char*) malloc((strlen(from->chrom)+1) * sizeof(char));
    strcpy( target->chrom, from->chrom);
    //copy the position level data
    target->pos=from->pos;
    target->ref_nuq=from->ref_nuq;
    target->n_samples=from->n_samples;
    ///copy the mutation related data
    target->mut_base=from->mut_base;
    target->mut_sample_idx=from->mut_sample_idx;
    target->mut_score=from->mut_score;
    target->mut_freq=from->mut_freq;
    target->cleanliness=from->cleanliness;
    strcpy( target->mut_type, from->mut_type);
    strncpy( target->mut_indel, from->mut_indel,MAX_INDEL_LEN);
    //loop over sample level data
    int i,j;
    for(i=0;i<target->n_samples;i++){
        //copy basic data cov, bases, quals
        target->raw_cov[i]=from->raw_cov[i];
        target->raw_bases[i] = (char*) malloc((strlen(from->raw_bases[i])+1) * sizeof(char));
        strcpy( target->raw_bases[i], from->raw_bases[i]);
        target->raw_quals[i] = (char*) malloc((strlen(from->raw_quals[i])+1) * sizeof(char));
        strcpy( target->raw_quals[i], from->raw_quals[i]);
        //copy filtered data, counts, coverage
        for(j=0;j<MAX_IDX;j++){
            target->counts[i][j]=from->counts[i][j];
            target->freqs[i][j]=from->freqs[i][j];
        }
        //copy all indel sequences
        target->ins_bases[i] = (char**) malloc(target->counts[i][INS_START_IDX] * sizeof(char*));
        target->del_bases[i] = (char**) malloc(target->counts[i][DEL_START_IDX] * sizeof(char*));
        for (j=0;j<target->counts[i][INS_START_IDX];j++){
            target->ins_bases[i][j] = (char*) malloc((strlen(from->ins_bases[i][j])+1) * sizeof(char));
            strcpy( target->ins_bases[i][j], from->ins_bases[i][j]);
        }
        for (j=0;j<target->counts[i][DEL_START_IDX];j++){
            target->del_bases[i][j] = (char*) malloc((strlen(from->del_bases[i][j])+1) * sizeof(char));
            strcpy( target->del_bases[i][j], from->del_bases[i][j]);
        }
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// formatted printing functions
////////////////////////////////////////////////////////////////////////////


/*
    printf formatted representation of mpileup line
*/
int print_mplp(struct mplp* my_mplp){
    //print position level info
    printf("%s %d %c\n",my_mplp->chrom,my_mplp->pos,my_mplp->ref_nuq);

    //print counts and freqs and bases and quals
    int i;
    for(i=0;i<my_mplp->n_samples;i++){
        printf("cov %d,A %d,C %d,G %d,T %d,del %d,ins_start %d, del_start %d,read_start %d, read_end %d\n",
               my_mplp->counts[i][COV_IDX],
               my_mplp->counts[i][A_IDX],
               my_mplp->counts[i][C_IDX],
               my_mplp->counts[i][G_IDX],
               my_mplp->counts[i][T_IDX],
               my_mplp->counts[i][DEL_IDX],
               my_mplp->counts[i][INS_START_IDX],
               my_mplp->counts[i][DEL_START_IDX],
               my_mplp->counts[i][READ_START_IDX],
               my_mplp->counts[i][READ_END_IDX]);
        printf("A %.2f,C %.2f,G %.2f,T %.2f,del %.2f,ins_start %.2f, del_start %.2f,read_start %.2f, read_end %.2f\n",
               my_mplp->freqs[i][A_IDX],
               my_mplp->freqs[i][C_IDX],
               my_mplp->freqs[i][G_IDX],
               my_mplp->freqs[i][T_IDX],
               my_mplp->freqs[i][DEL_IDX],
               my_mplp->freqs[i][INS_START_IDX],
               my_mplp->freqs[i][DEL_START_IDX],
               my_mplp->freqs[i][READ_START_IDX],
               my_mplp->freqs[i][READ_END_IDX]);
        printf("%d %s %s\n",my_mplp->raw_cov[i],my_mplp->raw_bases[i],my_mplp->raw_quals[i]);
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////
// process input mpileup line from samtools mpileup command output
////////////////////////////////////////////////////////////////////////////

/*
      process input mpileup line from samtools mpileup command output
*/
int process_mplp_input_line(struct mplp* my_mplp,char* line, ssize_t line_size,
                            int baseq_limit,char** sample_names,int n_sample_names){

    //read mpileup
    get_mplp(my_mplp,line,line_size);

    //assert if n_sampls==n_sample_names
    if(n_sample_names == my_mplp->n_samples){
        my_mplp->n_sample_names=n_sample_names;
        my_mplp->sample_names=sample_names;
    }
    else{
        printf("ERROR, length of sample names != number of samples in input stream\n ");
        exit(1);
    }

    //count bases
    count_bases_all_samples(my_mplp,baseq_limit);

    //calculate freqs
    calculate_freqs_all_samples(my_mplp);

    //collect indels
    collect_indels_all_samples(my_mplp,baseq_limit);

    return 0;
}



////////////////////////////////////////////////////////////////////////////
// read pileup struct from mpileup line
////////////////////////////////////////////////////////////////////////////

/*
    gets mpileup line from the char* line
*/
int get_mplp(struct mplp* my_mplp,char* line, ssize_t line_size){
    //store the raw line too
    free(my_mplp->raw_line);
    my_mplp->raw_line = (char*)malloc( (line_size+1) * sizeof(char));
    strcpy(my_mplp->raw_line,line);

    //temp buffer for reading those entries which will be formatted as not strings
    char* tmp_str=NULL;

    ssize_t i=0;
    while(i<line_size){
        //chrom
        get_next_entry(line,line_size,&i,&(my_mplp->chrom));
        //position
        get_next_entry(line,line_size,&i,&tmp_str);
        my_mplp->pos= (int) strtol(tmp_str,NULL,10);
        //ref nuq
        get_next_entry(line,line_size,&i,&tmp_str);
        my_mplp->ref_nuq=tmp_str[0];

        //read samples
        int temp_sample=0;
        while(i<line_size){
            //coverage
            get_next_entry(line,line_size,&i,&tmp_str);
            my_mplp->raw_cov[temp_sample]= (int) strtol(tmp_str,NULL,10);
            //bases
            get_next_entry(line,line_size,&i,&(my_mplp->raw_bases[temp_sample]));
            //quals
            get_next_entry(line,line_size,&i,&(my_mplp->raw_quals[temp_sample]));
            temp_sample++;
        }
        //store the number of samples
        my_mplp->n_samples=temp_sample;
    }
    //free resources
    free(tmp_str);
    return 0;
}

/*
    gets next entry from tab separated input string as null terminated char*
*/
int get_next_entry(char* line, ssize_t line_size, ssize_t* pointer, char** result){
    int c,c0;
    c0 = (int) *pointer;
    while(line[*pointer]!='\t' && *pointer<line_size)(*pointer)++;
    c = (int) *pointer-c0;
    (*pointer)++;

    if(*result != NULL) free(*result);
    *result = (char*)malloc( (c+1) * sizeof(char));
    memcpy(*result,line+c0,c * sizeof(char));
    (*result)[c] = 0;

    return c;
}



////////////////////////////////////////////////////////////////////////////
// Count bases in mplieup struct
////////////////////////////////////////////////////////////////////////////

/*
    counts bases,indels, and read start end signs in all samples
*/
int count_bases_all_samples(struct mplp* my_mplp,int baseq_lim){
    int i=0;
    for(i=0;i<my_mplp->n_samples;i++){
        count_bases(my_mplp->raw_bases[i],my_mplp->raw_quals[i],my_mplp->counts[i],
                    my_mplp->ref_nuq,baseq_lim);
    }
    return 0;
}

/*
    counts bases in one sample
*/
int count_bases(char* bases, char* quals,int* counts,char ref_base,int baseq_lim){
    //initialize counts to zero
    int i,j;
    for( i=0;i<MAX_IDX;i++) counts[i]=0;

    i = 0; //pointer in str for bases
    j = 0; //pointer in str for qualities
    counts[COV_IDX]=0;
    while(bases[i]!=0){
        //beginning and end of the read signs
        if(bases[i] == '^' ) {i+=2;counts[READ_START_IDX]++;}
        else if(bases[i] == '$' ) {i++,counts[READ_END_IDX]++;}
        //deletions
        else if(bases[i]=='-' ) handle_deletion(bases,&counts[DEL_START_IDX],&i,quals[j-1],baseq_lim);
        //insetions
        else if(bases[i]=='+' ) handle_insertion(bases,&counts[INS_START_IDX],&i,quals[j-1],baseq_lim);
        //real base data
        else handle_base(bases,quals,counts,&counts[COV_IDX],&i,&j,baseq_lim);
    }
    //add refbase to corresponding base
    if(ref_base=='A') counts[A_IDX]+=counts[REF_IDX];
    else if(ref_base=='C') counts[C_IDX]+=counts[REF_IDX];
    else if(ref_base=='G') counts[G_IDX]+=counts[REF_IDX];
    else if(ref_base=='T') counts[T_IDX]+=counts[REF_IDX];
    return 0;
}


/*
    parse a base from the bases and quals
*/
int handle_base(char* bases,char* quals,int* base_counts, int* filtered_cov,
                   int* base_ptr,int* qual_ptr,int baseq_lim){

    char c = bases[*base_ptr];
    if(quals[*qual_ptr] >= baseq_lim + 33 ){
        if(c=='.' || c==',' )      base_counts[REF_IDX]++;
        else if(c=='A' || c=='a' ) base_counts[A_IDX]++;
        else if(c=='C' || c=='c' ) base_counts[C_IDX]++;
        else if(c=='G' || c=='g' ) base_counts[G_IDX]++;
        else if(c=='T' || c=='t' ) base_counts[T_IDX]++;
        else if(c=='*' ) base_counts[DEL_IDX]++;
        (*filtered_cov)++;
    }
    (*qual_ptr)++;
    (*base_ptr)++;

    return 0;
}

/*
    parse a deletion from the bases and quals
*/
int handle_deletion(char* bases,int* del_count,int* base_ptr,char qual,int baseq_lim){
    char* offset;
    int indel_len= (int) strtol(&bases[*base_ptr+1],&offset,10);
    (*base_ptr)+= offset-&bases[*base_ptr] + indel_len;
    if(qual >= baseq_lim + 33 )(*del_count)++;
    return 0;
}

/*
    parse an insertion from the bases and quals
*/
int handle_insertion(char* bases,int* ins_count,int* base_ptr,char qual,int baseq_lim){
    char* offset;
    int indel_len= (int) strtol(&bases[*base_ptr+1],&offset,10);
    (*base_ptr)+= offset-&bases[*base_ptr] + indel_len;
    if(qual >= baseq_lim + 33 ) (*ins_count)++;
    return 0;
}


///////////////////////////////////////////////////////////////////////////
// Calculate base frequencies
////////////////////////////////////////////////////////////////////////////

/*
    calculate base freqs in all samples
*/
int calculate_freqs_all_samples(struct mplp* my_mplp){
    int i=0;
    for(i=0;i<my_mplp->n_samples;i++){
        calculate_freqs(my_mplp->freqs[i],my_mplp->counts[i]);
    }
    return 0;
}


/*
   calculate freqs in a sample
*/
int calculate_freqs(double* freqs,int* counts){
    int i;
    if (counts[COV_IDX]!=0){
        for(i=0;i<MAX_IDX;i++){
            freqs[i]=( (double) counts[i]) / counts[COV_IDX];
        }
    }
    else {
        for(i=0;i<MAX_IDX;i++){
            freqs[i]= ZERO_COV_FREQ;
        }
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// Collect indels
////////////////////////////////////////////////////////////////////////////

/*
    collect indels in all samples
*/
int collect_indels_all_samples(struct mplp* my_mplp,int baseq_lim){
    int i;
    for(i=0;i<my_mplp->n_samples;i++){
        collect_indels(my_mplp->raw_bases[i],
                       my_mplp->raw_quals[i],
                       &my_mplp->ins_bases[i],
                       my_mplp->counts[i][INS_START_IDX],
                       &my_mplp->del_bases[i],
                       my_mplp->counts[i][DEL_START_IDX],
                       baseq_lim);
    }
    return 0;
}

/*
    collect the inserted, and deleted bases
*/
int collect_indels(char* bases,char* quals, char*** ins_bases, int ins_count,
                   char*** del_bases,int del_count, int baseq_lim){
    //allocate new memory
    if( *ins_bases!=NULL || *del_bases!=NULL){
        printf("ERROR: collect_indels() called on not NULL ins_bases, del_bases pointers,\n");
        printf("       maybe its called 2nd time, or mplp struct not freed, or not initialized,\n");
        printf("       possible memory leak, exiting");
        exit(1);
    }
    *ins_bases = (char**) malloc( (ins_count) * sizeof(char*));
    *del_bases = (char**) malloc( (del_count) * sizeof(char*));

    int i,j,del_c,ins_c; //pointers in data
    i = j = del_c = ins_c = 0;
    char* offset;
    while(bases[i]!=0){
        //beginning and end of the read signs
        if(bases[i] == '$' ) i++; //next
        else if(bases[i] == '^' ) i+=2; //jump next character (mapq too)
        //deletions
        else if(bases[i]=='-' ) {
            int indel_len= (int) strtol(&bases[i+1],&offset,10);
            i+= offset-&bases[i] + indel_len;
            if( quals[j-1] >= baseq_lim + 33 ){
                (*del_bases)[del_c] = (char*) malloc( (indel_len+1) * sizeof(char));
                memcpy((*del_bases)[del_c],offset,indel_len * sizeof(char));
                (*del_bases)[del_c][indel_len]=0;
                del_c++;
            }
        }
        //insertions
        else if(bases[i]=='+' ){
            int indel_len= (int) strtol(&bases[i+1],&offset,10);
            i+= offset-&bases[i] + indel_len;
            if( quals[j-1] >= baseq_lim + 33 ){
                (*ins_bases)[ins_c] = (char*) malloc( (indel_len+1) * sizeof(char));
                memcpy((*ins_bases)[ins_c],offset,indel_len * sizeof(char));
                (*ins_bases)[ins_c][indel_len]=0;
                ins_c++;
            }
        }
        //real base data
        else {i++;j++;}
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////
// Print 
////////////////////////////////////////////////////////////////////////////

/*
    print base freqs
*/
int print_freq(struct mplp* my_mplp){
    int i=0;
    for(i=0;i<my_mplp->n_samples-1;i++){
        //printf("%.2f\t",my_mplp->freqs[i][REF_IDX]);
        printf("%.2f\t",my_mplp->freqs[i][A_IDX]);
        printf("%.2f\t",my_mplp->freqs[i][C_IDX]);
        printf("%.2f\t",my_mplp->freqs[i][G_IDX]);
        printf("%.2f\t",my_mplp->freqs[i][T_IDX]);
    }
    //printf("%.2f\t",my_mplp->freqs[i][REF_IDX]);
    printf("%.2f\t",my_mplp->freqs[i][A_IDX]);
    printf("%.2f\t",my_mplp->freqs[i][C_IDX]);
    printf("%.2f\t",my_mplp->freqs[i][G_IDX]);
    printf("%.2f\n",my_mplp->freqs[i][T_IDX]);
    return 0;
}

/*
 print base counts
 */
int print_counts(struct mplp* my_mplp){
    int i=0;
    for(i=0;i<my_mplp->n_samples-1;i++){
        //printf("%d\t",my_mplp->counts[i][REF_IDX]);
        printf("%d\t",my_mplp->counts[i][A_IDX]);
        printf("%d\t",my_mplp->counts[i][C_IDX]);
        printf("%d\t",my_mplp->counts[i][G_IDX]);
        printf("%d\t",my_mplp->counts[i][T_IDX]);
    }
    //printf("%d\t",my_mplp->counts[i][REF_IDX]);
    printf("%d\t",my_mplp->counts[i][A_IDX]);
    printf("%d\t",my_mplp->counts[i][C_IDX]);
    printf("%d\t",my_mplp->counts[i][G_IDX]);
    printf("%d\n",my_mplp->counts[i][T_IDX]);
    return 0;
}
