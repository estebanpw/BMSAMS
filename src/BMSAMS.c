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
    char seq_A[MAX_READ_SIZE] = "CTTAGATCGTACCAAAATATTAC";
    char seq_B[MAX_READ_SIZE] = "ATTAGATCGTACCACTATAAGTTTAC";
    char seq_C[MAX_READ_SIZE] = "CTTAGATCGTTCCACACATATAC";
    char seq_D[MAX_READ_SIZE] = "CTTAGATCGTACCACACAATTAC";

    char reco_P1_A[MAX_READ_SIZE]; 
    char reco_P1_B[MAX_READ_SIZE]; 
    char reco_P2_A[MAX_READ_SIZE]; 
    char reco_P2_B[MAX_READ_SIZE]; 
    char reco_P3_A[MAX_READ_SIZE];
    char reco_P3_B[MAX_READ_SIZE];
    


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
    
    uint64_t xlen;
    uint64_t ylen;

    /*void build_multiple_alignment(char * reconstruct_X, char * reconstruct_Y, char * my_x, char * my_y, struct cell ** table, struct positioned_cell * mc, char * writing_buffer_alignment, uint64_t xlen, uint64_t ylen, int64_t * cell_path_y, long double * window, int64_t iGap, int64_t eGap, BasicAlignment * ba)*/
    
    BasicAlignment ba;
    seq_piles_up_X[0] = seq_A;
    seq_piles_up_Y[0] = seq_B;
    build_unanchored_alignment(cell_path_y, seq_piles_up_X, seq_piles_up_Y, 1, 1);
    uint64_t xlen = get_max_length_of_sequences(seq_piles_up_X, 1);
    uint64_t ylen = get_max_length_of_sequences(seq_piles_up_Y, 1);
    Two_seqs ts = build_multiple_alignment(reco_P1_A, reco_P1_B, seq_piles_up_X, seq_piles_up_Y, 1, 1, table, mc, writing_buffer_alignment, xlen, ylen, cell_path_y, &window, iGap, eGap, &ba);

    fprintf(stdout, "%s\n", ts.s1);
    fprintf(stdout, "%s\n", ts.s2);

    /*

    fprintf(stdout, " ---- \n");
    build_unanchored_alignment(cell_path_y, strlen(seq_C), strlen(seq_D));
    ts = build_multiple_alignment(reco_P2_A, reco_P2_B, seq_C, seq_D, table, mc, writing_buffer_alignment, strlen(seq_C), strlen(seq_D), cell_path_y, &window, iGap, eGap, &ba);

    fprintf(stdout, "%s\n", ts.s1);
    fprintf(stdout, "%s\n", ts.s2);

    fprintf(stdout, " ---- \n");
    //build_unanchored_alignment(cell_path_y, MAX(strlen(reco_P1_A), strlen(reco_P1_B)), MAX(strlen(reco_P2_A), strlen(reco_P2_B)));
    //ts = build_multiple_alignment(reco_P3_A, reco_P3_B, reco_P2_A, seq_D, table, mc, writing_buffer_alignment, strlen(seq_C), strlen(seq_D), cell_path_y, &window, iGap, eGap, &ba);
    */
    return 0;    
}



