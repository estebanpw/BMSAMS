/*********

File        BMSAMS.c
Author      EPW <estebanpw@uma.es>
Description This is just a file to call the algorithm. Use it as a pipelined program.


**********/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) <= (y)) ? (x) : (y))
#define STARTING_SEQS 1000
#define PIECE_OF_DB_REALLOC 3200000 //half a gigabyte if divided by 8 bytes



int main(int argc, char ** av){

    uint64_t i;
    int64_t iGap = -3, eGap = -1;
    /*
    // Example seqs 1
    char seq_A[MAX_READ_SIZE] = "CTTAGATCGTACCAAAATATTACCCCCCTCGT";
    char seq_B[MAX_READ_SIZE] = "ATTAGATCGTACCACTATAAGTTTACCCCCCAAAATCGT";
    char seq_C[MAX_READ_SIZE] = "CTTAGATCGTTCCACACATATACCCCCCACTCGT";
    char seq_D[MAX_READ_SIZE] = "CTTAGATCGTACCACACAATTACCCCCCACTCGT";
    */
    // Example seqs 2
    char seq_A[MAX_READ_SIZE] = "CTTAGATCGTACCAAAATTTTTTATTACCCCCCTCGT";
    char seq_B[MAX_READ_SIZE] = "CTTAGATCGTACCAAAATATTACCCCCCTCGT";
    char seq_C[MAX_READ_SIZE] = "CTTAGATCGTACCAAAATATTACCCCCCTCGT";
    char seq_D[MAX_READ_SIZE] = "CTTAGATCGTACCAATTTTTAATATTACCCCCCTCGT";

    printf("Init seqs\n");
    printf("%s\n", seq_A);
    printf("%s\n", seq_C);
    printf("%s\n", seq_B);
    printf("%s\n-------------\n", seq_D);

    char reco_P1_A[MAX_READ_SIZE]; 
    char reco_P1_B[MAX_READ_SIZE]; 
    char reco_P2_A[MAX_READ_SIZE]; 
    char reco_P2_B[MAX_READ_SIZE]; 
    char reco_P3_A[MAX_READ_SIZE];
    char reco_P3_B[MAX_READ_SIZE];
    char aux[MAX_READ_SIZE];
    


    struct positioned_cell * mc = (struct positioned_cell *) malloc(MAX_READ_SIZE * sizeof(struct positioned_cell));
    struct cell ** table = (struct cell **) malloc(MAX_READ_SIZE * sizeof(struct cell *));
    char * writing_buffer_alignment = (char *) malloc(2*MAX_READ_SIZE * sizeof(char));
    uint64_t j;
    for(j=0;j<MAX_READ_SIZE;j++){
        table[j] = (struct cell *) malloc(MAX_WINDOW_SIZE*sizeof(struct cell));
    }
    long double window = 0.5;

    // Generalization part here 

    int64_t * cell_path_y = (int64_t *) malloc(MAX_READ_SIZE*sizeof(int64_t));
    if(cell_path_y == NULL) terror("Could not allocate cell path");



    char ** seq_piles_up_X = (char **) malloc(10 * sizeof(char *));
    char ** seq_piles_up_Y = (char **) malloc(10 * sizeof(char *));

    char ** recons_X = (char **) malloc(10 * sizeof(char *));
    char ** recons_Y = (char **) malloc(10 * sizeof(char *));
    

    /*void build_multiple_alignment(char * reconstruct_X, char * reconstruct_Y, char * my_x, char * my_y, struct cell ** table, struct positioned_cell * mc, char * writing_buffer_alignment, uint64_t xlen, uint64_t ylen, int64_t * cell_path_y, long double * window, int64_t iGap, int64_t eGap, BasicAlignment * ba)*/
    
    BasicAlignment ba;
    
    seq_piles_up_X[0] = seq_A;
    seq_piles_up_Y[0] = seq_B;
    recons_X[0] = reco_P1_A;
    recons_Y[0] = reco_P1_B;



    build_unanchored_alignment(cell_path_y, seq_piles_up_X, seq_piles_up_Y, 1, 1);
    uint64_t xlen = get_max_length_of_sequences(seq_piles_up_X, 1);
    uint64_t ylen = get_max_length_of_sequences(seq_piles_up_Y, 1);
    All_seqs as = build_multiple_alignment(recons_X, recons_Y, seq_piles_up_X, seq_piles_up_Y, 1, 1, table, mc, writing_buffer_alignment, xlen, ylen, cell_path_y, &window, iGap, eGap, &ba, aux);

    printf("alignment so far\n");
    printf("X 0 -> %s (%"PRIu64")\n", seq_piles_up_X[0], strlen(seq_piles_up_X[0]));
    printf("Y 0 -> %s (%"PRIu64")\n", seq_piles_up_Y[0], strlen(seq_piles_up_Y[0]));
        

    fprintf(stdout, " ---- \n");
    
    seq_piles_up_X[1] = seq_C;
    seq_piles_up_Y[1] = seq_D;
    recons_X[1] = reco_P2_A;
    recons_Y[1] = reco_P2_B;

    char ** other_alignment_X = &seq_piles_up_X[1];
    char ** other_alignment_Y = &seq_piles_up_Y[1];
    char ** other_recons_X = &recons_X[1];
    char ** other_recons_Y = &recons_Y[1];

    build_unanchored_alignment(cell_path_y, other_alignment_X, other_alignment_Y, 1, 1);
    xlen = get_max_length_of_sequences(other_alignment_X, 1);
    ylen = get_max_length_of_sequences(other_alignment_Y, 1);
    as = build_multiple_alignment(other_recons_X, other_recons_Y, other_alignment_X, other_alignment_Y, 1, 1, table, mc, writing_buffer_alignment, xlen, ylen, cell_path_y, &window, iGap, eGap, &ba, aux);
    
    printf("alignment so far\n");
    printf("X 1 -> %s\n", seq_piles_up_X[1]);
    printf("Y 1 -> %s\n", seq_piles_up_Y[1]);

    fprintf(stdout, " ---- \n");

    char ** group_X = (char **) malloc(10 * sizeof(char *));
    char ** group_Y = (char **) malloc(10 * sizeof(char *));
    char ** g_reco_X = (char **) malloc(10 * sizeof(char *));
    char ** g_reco_Y = (char **) malloc(10 * sizeof(char *));

    group_X[0] = seq_piles_up_X[0];
    group_X[1] = seq_piles_up_Y[0];
    group_Y[0] = seq_piles_up_X[1];
    group_Y[1] = seq_piles_up_Y[1];

    g_reco_X[0] = recons_X[0];
    g_reco_X[1] = recons_Y[0];
    g_reco_Y[0] = recons_X[1];
    g_reco_Y[1] = recons_Y[1];

    build_unanchored_alignment(cell_path_y, group_X, group_Y, 2, 2);
    xlen = get_max_length_of_sequences(group_X, 2);
    ylen = get_max_length_of_sequences(group_Y, 2);
    as = build_multiple_alignment(g_reco_X, g_reco_Y, group_X, group_Y, 2, 2, table, mc, writing_buffer_alignment, xlen, ylen, cell_path_y, &window, iGap, eGap, &ba, aux);

    printf("Lengths final %"PRIu64", %"PRIu64"\n", xlen, ylen);

    printf("Final alignment\n");
    printf("X 0 -> %s\n", seq_piles_up_X[0]);
    printf("Y 0 -> %s\n", seq_piles_up_Y[0]);
    printf("X 1 -> %s\n", seq_piles_up_X[1]);
    printf("Y 1 -> %s\n", seq_piles_up_Y[1]);

    //build_unanchored_alignment(cell_path_y, MAX(strlen(reco_P1_A), strlen(reco_P1_B)), MAX(strlen(reco_P2_A), strlen(reco_P2_B)));
    //ts = build_multiple_alignment(reco_P3_A, reco_P3_B, reco_P2_A, seq_D, table, mc, writing_buffer_alignment, strlen(seq_C), strlen(seq_D), cell_path_y, &window, iGap, eGap, &ba);
    
    return 0;    
}



