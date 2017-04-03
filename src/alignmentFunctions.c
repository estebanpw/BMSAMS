#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>
#include "structs.h"
#include "alignmentFunctions.h"
#include "commonFunctions.h"

void pile_chars_up(char ** list_seq_X, char ** list_seq_Y, char * list_piled_X, char * list_piled_Y, uint64_t nx, uint64_t ny, uint64_t posx, uint64_t posy){
    uint64_t i;
    for(i=0;i<nx;i++) list_piled_X[i] = list_seq_X[i][posx];
    for(i=0;i<ny;i++) list_piled_Y[i] = list_seq_Y[i][posy];
}

int64_t compare_piled_up_chars(char * list_piled_X, char * list_piled_Y, uint64_t nx, uint64_t ny){

    if(nx == 1 && ny == 1) return compare_letters(list_piled_X[0], list_piled_Y[0]);

    uint64_t i;
    int64_t acum = 0;
    uint64_t types_X[6] = {0,0,0,0,0,0}; // A, C, G, T, N, -
    uint64_t types_Y[6] = {0,0,0,0,0,0};
    for(i=0;i<nx;i++){
        switch(list_piled_X[i]){
            case 'A': types_X[0]++; break;
            case 'C': types_X[1]++; break;
            case 'G': types_X[2]++; break;
            case 'T': types_X[3]++; break;
            case 'N': types_X[4]++; break;
            case '-': types_X[5]++; break;
        }
    }
    for(i=0;i<ny;i++){
        switch(list_piled_Y[i]){
            case 'A': types_Y[0]++; break;
            case 'C': types_Y[1]++; break;
            case 'G': types_Y[2]++; break;
            case 'T': types_Y[3]++; break;
            case 'N': types_Y[4]++; break;
            case '-': types_Y[5]++; break;
        }
    }
    for(i=0;i<4;i++){
        if(types_X[i] > 0 && types_Y[i] > 0) acum += (types_X[i] + types_Y[i])*POINT;
        types_X[i] -= types_Y[i];
        acum += labs(types_X[i]) * (-POINT);
    }
    for(i=4;i<6;i++){
        acum += (types_X[i] + types_Y[i])*(-POINT);
    }

    return acum;
}

int64_t compare_letters(char a, char b){
    if(a != 'N' && a != '-') return (a == b) ? POINT : -POINT;
    return -POINT;
}

uint64_t get_max_length_of_sequences(char ** a, uint64_t n_x){
    uint64_t max_len = 0;
    uint64_t len;
    uint64_t i;
    for(i=0;i<n_x;i++){
        len = strlen(a[i]);
        if(len > max_len) max_len = len;
    }
    return max_len;
}


void build_unanchored_alignment(int64_t * cell_path_y, char ** seq_piles_up_X, char ** seq_piles_up_Y, uint64_t n_x, uint64_t n_y){

    uint64_t l_A = get_max_length_of_sequences(seq_piles_up_X, n_x), l_B = get_max_length_of_sequences(seq_piles_up_Y, n_y);
    
    Point p0 = {0,0};
    Point p1 = {l_A/4, l_B/4};
    Point p2 = {l_A - l_A/4, l_B - l_B/4};
    Point p3 = {l_A, l_B};
    calculate_y_cell_path(p0, p1, p2, p3, cell_path_y);
}



All_seqs build_multiple_alignment(char ** reconstruct_X, char ** reconstruct_Y, char ** my_x, char ** my_y, uint64_t nx, uint64_t ny, struct cell ** table, struct positioned_cell * mc, char * writing_buffer_alignment, uint64_t xlen, uint64_t ylen, int64_t * cell_path_y, long double * window, int64_t iGap, int64_t eGap, BasicAlignment * ba, char * aux){
 

    //Do some printing of alignments here
    uint64_t i, j, k, curr_window_size;
    All_seqs as;
    
    struct best_cell bc = NW(my_x, 0, xlen, my_y, 0, ylen, iGap, eGap, table, mc, 0, cell_path_y, window, &curr_window_size, nx, ny);

    uint64_t max_len = MAX(xlen, ylen);

    for(k=0;k<nx;k++){
        backtrackingNW(my_x[k], 0, xlen, my_x[k], 0, ylen, table, reconstruct_X[k], aux, &bc, &i, &j, cell_path_y, curr_window_size, ba);
        i++; j++;
        fprintf(stdout, "%s\n", &reconstruct_X[k][i]);
        memcpy(&my_x[k][0], &reconstruct_X[k][i], max_len);
    }


    for(k=0;k<ny;k++){
        backtrackingNW(my_y[k], 0, xlen, my_y[k], 0, ylen, table, aux, reconstruct_Y[k], &bc, &i, &j, cell_path_y, curr_window_size, ba);
        i++; j++;    
        fprintf(stdout, "%s\n", &reconstruct_Y[k][j]);
        memcpy(&my_y[k][0], &reconstruct_Y[k][j], max_len);
    }

    
    //backtrackingNW(my_x, 0, xlen, my_y, 0, ylen, table, reconstruct_X, reconstruct_Y, &bc, &i, &j, cell_path_y, curr_window_size, ba, nx, ny);
    //uint64_t offset = 0, before_i = 0, before_j = 0;
    //i++; j++;

    


    

    

    /*
    Two_seqs ts;
    ts.s1 = &reconstruct_X[i];
    ts.s2 = &reconstruct_Y[j];
    */

    return as;

}

void calculate_y_cell_path(Point p0, Point p1, Point p2, Point p3, int64_t * y_points){
    
    //Calculate lines between points
    uint64_t i;

    #ifdef VERBOSE
    printf("Built on\n");
    printf("(%"PRIu64", %"PRIu64")\n", p0.x, p0.y);
    printf("(%"PRIu64", %"PRIu64")\n", p1.x, p1.y);
    printf("(%"PRIu64", %"PRIu64")\n", p2.x, p2.y);
    printf("(%"PRIu64", %"PRIu64")\n", p3.x, p3.y);
    #endif

    long double deltax, deltay, deltaerr, error;
    uint64_t y;

    //P0 to P1
    deltax = p1.x - p0.x;
    deltay = p1.y - p0.y;
    if(deltax != 0) deltaerr = fabsl(deltay/deltax); else deltaerr = 0;
    //printf("Deltas  x: %Le y: %Le Error: %Le\n", deltax, deltay, deltaerr);
    error = deltaerr - 0.5;
    y = p0.y;

    for(i=p0.x;i<p1.x;i++){
        y_points[i] = (int64_t) y;
        error = error + deltaerr;
        if(error >= 0.5){
            y++;
            error = error - 1;
        }
    }

    //P1 to P2

    deltax = p2.x - p1.x;
    deltay = p2.y - p1.y;
    if(deltax != 0) deltaerr = fabsl(deltay/deltax); else deltaerr = 0;
    //printf("Deltas  x: %Le y: %Le Error: %Le\n", deltax, deltay, deltaerr);
    error = deltaerr - 0.5;
    y = p1.y;

    for(i=p1.x;i<p2.x;i++){
        y_points[i] = (int64_t) y;
        error = error + deltaerr;
        if(error >= 0.5){
            y++;
            error = error - 1;
        }
    }
    
    //P2 to P3

    deltax = p3.x - p2.x;
    deltay = p3.y - p2.y;
    if(deltax != 0) deltaerr = fabsl(deltay/deltax); else deltaerr = 0;
    //printf("Deltas  x: %Le y: %Le Error: %Le\n", deltax, deltay, deltaerr);
    error = deltaerr - 0.5;
    y = p2.y;

    for(i=p2.x;i<p3.x;i++){
        y_points[i] = (int64_t) y;
        error = error + deltaerr;
        if(error >= 0.5){
            y++;
            error = error - 1;
        }
    }

    #ifdef VERBOSE
    for(i=0;i<p3.x;i++){
        printf("%"PRIu64" -> ", y_points[i]);
        if(i % 50 == 0) getchar();
    }
    #endif
    

}

struct best_cell NW(char ** X, uint64_t Xstart, uint64_t Xend, char ** Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell ** table, struct positioned_cell * mc, int show, int64_t * cell_path_y, long double * window, uint64_t * current_window_size, uint64_t nx, uint64_t ny){
    

    uint64_t i, j, j_prime;
    int64_t scoreDiagonal,scoreLeft,scoreRight,score, delta_dif_1_0, delta_dif_2_1, limit_left, limit_right, j_right_prime = 1, j_left_prime, j_diag_prime;

    struct best_cell bc;
    bc.c.score = INT64_MIN;

    // To pile up chars 
    char list_piled_X[nx];
    char list_piled_Y[ny];
    
    //The window size will be a +-15% of the square root of the product of lengths
    int64_t window_size = (uint64_t) (*window * sqrtl((long double) Xend * (long double) Yend));
    *current_window_size = (uint64_t) window_size;

    //The limits to the window
    limit_left = 0;
    limit_right = 2*window_size + 1;
    if(limit_right > MAX_WINDOW_SIZE) limit_right = MAX_WINDOW_SIZE;
    
    struct positioned_cell mf;
    mf.score = INT64_MIN;
    

    //First row. iCounter serves as counter from zero
    //printf("..0%%");
    //Zero will always be
    pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, 0, 0);
    table[0][0].score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny);
    //table[0][0].score = compare_letters(X[0], Y[0]);
    mc[0].score = table[0][0].score;
    mc[0].xpos = 0;
    mc[0].ypos = 0;

    for(i=1;i<Yend;i++){
        //table[0][i].score = (X[0] == Y[i]) ? POINT : -POINT;
        if(i < cell_path_y[0] + window_size){
            pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, 0, i);
            table[0][i].score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny) + iGap + (i-1)*eGap;
            //table[0][i].score = compare_letters(X[0], Y[i]) + iGap + (i-1)*eGap;
        } 
        //table[Xstart][i].xfrom = Xstart;
        //table[Xstart][i].yfrom = i;
        //Set every column max
        
        mc[i].score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny) + iGap + (i-1)*eGap;
        #ifdef VERBOSE
        printf("%02"PRId64" ", mc[i].score);
        #endif
        mc[i].xpos = 0;
        mc[i].ypos = i;

    }
    #ifdef VERBOSE
    printf("\n");
    #endif
    //Set row max
    mf.score = table[0][0].score;
    mf.xpos = 0;
    mf.ypos = 0;
    //Init j
    j = MAX(1,(cell_path_y[1] - window_size));

    //Go through full matrix
    for(i=1;i<Xend;i++){
        //Fill first rowcell

        //Conversion for the j-coordinate
        j_prime = 1;

        //table[i][0].score = (X[i] == Y[0]) ? POINT : -POINT;
        if(cell_path_y[i] - window_size <= 0){
            pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, i, 0);
            table[i][0].score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny) + iGap + (i-1)*eGap;
            //table[i][0].score = compare_letters(X[i], Y[0]) + iGap + (i-1)*eGap;
            mf.score = table[i][0].score;
        }else{
            pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, i, 0);
            mf.score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny) + iGap + (i-1)*eGap;
        }

        mf.xpos = i-1;
        mf.ypos = 0;

        delta_dif_1_0 = MAX(1, (cell_path_y[i] - window_size)) - MAX(1,(cell_path_y[i-1] - window_size)); //j-1
        if(i>1) delta_dif_2_1 = MAX(1, (cell_path_y[i-1] - window_size)) - MAX(1, (cell_path_y[i-2] - window_size)); //j-2

        #ifdef VERBOSE 
        printf("D1_0: %"PRId64" D2_1: %"PRId64"\n", delta_dif_1_0, delta_dif_2_1);
        #endif

        #ifdef VERBOSE
        printf("%02"PRId64" ", mf.score);
        #endif
        //printf("Check on i: (%"PRIu64") from - to (%"PRIu64", %"PRIu64")\n", i, 0L, Xend);
        //printf("I will go from %"PRIu64" to %"PRIu64"\n", (uint64_t) MAX(1,(cell_path_y[i] - window_size)), (uint64_t) MIN(Yend,(cell_path_y[i] + window_size)));
        //getchar();

        #ifdef VERBOSE
        int64_t r;
        for(r=0;r<MAX(0,(cell_path_y[i] - window_size)); r++){
            printf("  ");
        }
        #endif

        

        

        for(j=MAX(1,(cell_path_y[i] - window_size));j<MIN(Yend,(cell_path_y[i] + window_size));j++){
            //printf("Doing on : (%"PRIu64",%"PRIu64" and jprime=%"PRIu64"\n", i,j,j_prime);
            //Check if max in row has changed
            //if(j > MAX(1, cell_path_y[i-1] - window_size +1) && mf.score <= table[i][j-2].score){

            //Calculate the real j position in the windowed table
            
            j_left_prime = ((int64_t)j_prime - (2 - delta_dif_1_0));
            //j_diag_prime = ((int64_t)j_prime - (1 - delta_dif_1_0));
            j_diag_prime = ((int64_t)j_prime - (1 - delta_dif_1_0));
            if(i > 1){
                j_right_prime = ((int64_t)j_prime - (1 - (delta_dif_1_0 + delta_dif_2_1)));
            }

            if(j > MAX(1, cell_path_y[i-1] - window_size +1) && j_left_prime >= limit_left && j_left_prime < limit_right && table[i-1][j_left_prime].score >= mf.score){
                //mf.score = table[i-1][j-2].score;
                mf.score = table[i-1][j_left_prime].score;
                mf.xpos = i-1;
                mf.ypos = j-2;
            }
            //printf("RowMax: %"PRId64"@(%"PRIu64", %"PRIu64")\t", mf.score, mf.xpos, mf.ypos);
            
            //score = (X[i] == Y[j]) ? POINT : -POINT;
            pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, i, j);
            score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny);
            //score = compare_letters(X[i], Y[j]);

            //Precondition: Upper row needs to reach up to diagonal
            //if((cell_path_y[i-1]+window_size) >= j-1){
            if(j-1 >= MAX(0, (cell_path_y[i-1]-window_size)) && (cell_path_y[i-1]+window_size) >= j-1 && j_diag_prime >= limit_left && j_diag_prime < limit_right){
                //scoreDiagonal = table[i-1][j-1].score + score;
                //printf("prevdiag: %"PRId64"\n", table[i-1][j_diag_prime].score);
                scoreDiagonal = table[i-1][j_diag_prime].score + score;
                //printf("j_diag: %"PRId64":", j_diag_prime);
            }else{
                scoreDiagonal = INT64_MIN;
            }
            
            if(i>=1 && j>1){
                scoreLeft = mf.score + iGap + (j - (mf.ypos+2))*eGap + score;
                
            }else{
                scoreLeft = INT64_MIN;
            }

            if(j>=1 && i>1){
                scoreRight = mc[j-1].score + iGap + (i - (mc[j-1].xpos+2))*eGap + score;
                //if(scoreRight == -12) printf("MC: %"PRId64", From: %"PRIu64", %"PRIu64"->", mc[j-1].score, mc[j-1].xpos, mc[j-1].ypos);
            }else{
                scoreRight = INT64_MIN;
            }
            
            //Choose maximum
            /*
            #ifdef VERBOSE
            printf("The game starts at %"PRId64"\n", MAX(0, cell_path_y[i] - window_size));
            printf("from %c %c and I get to %"PRIu64" while j=%"PRIu64"\n", X[i], Y[j], j_prime, j);
            printf("j_prime: %"PRId64"\n", j_prime);
            printf("j_diag_prime: %"PRId64" limits[%"PRId64", %"PRId64"]\n", j_diag_prime, limit_left, limit_right);
            printf("Score DIAG: %"PRId64"; LEFT: %"PRId64"; RIGHT: %"PRId64"\n", scoreDiagonal, scoreLeft, scoreRight);
            printf("currmf: %"PRId64" mc: %"PRId64"\n", mf.score, mc[j-1].score);
            #endif
            */
            
            if(scoreDiagonal >= scoreLeft && scoreDiagonal >= scoreRight){
                //Diagonal
                table[i][j_prime].score = scoreDiagonal;
                table[i][j_prime].xfrom = i-1;
                table[i][j_prime].yfrom = j-1;
                
                                
            }else if(scoreRight > scoreLeft){
                table[i][j_prime].score = scoreRight;
                table[i][j_prime].xfrom = mc[j-1].xpos;
                table[i][j_prime].yfrom = mc[j-1].ypos;
                
            }else{
                table[i][j_prime].score = scoreLeft;
                table[i][j_prime].xfrom = mf.xpos;
                table[i][j_prime].yfrom = mf.ypos;
            }
        
        
            //check if column max has changed
            //New condition: check if you filled i-2, j-1
            
            if(i > 1 && j >= 1 && j_right_prime >= limit_left && j_right_prime < limit_right && table[i-2][j_right_prime].score >= mc[j-1].score){
                //mc[j-1].score = table[i-2][j-(1+j_prime)].score;
                //Should be the j_prime we had at cell_path_y
                
                mc[j-1].score = table[i-2][j_right_prime].score;
                mc[j-1].xpos = i-2;
                mc[j-1].ypos = j-1;
            }
            if(i == Xend-1 || j == Yend-1){

                if(i == Xend-1 && j != Yend-1){
            		table[i][j_prime].score = table[i][j_prime].score + iGap + (Yend - j)*eGap;
            	}else if(j == Yend-1 && i != Xend-1){
            		table[i][j_prime].score = table[i][j_prime].score + iGap + (Xend - i)*eGap;
            	}
                //Check for best cell
                if(table[i][j_prime].score >= bc.c.score){ bc.c.score = table[i][j_prime].score; bc.c.xpos = i; bc.c.ypos = j; bc.j_prime = j_prime; }
            }
            #ifdef VERBOSE
            //printf("Put score: %"PRId64"\n\n", table[i][j_prime].score);
            printf("(%"PRId64")%02"PRId64" ", j_diag_prime, table[i][j_prime].score); //printf("->(%"PRIu64", %"PRIu64")", i, j); printf("[%c %c]", X[i], Y[j]);
            //if(scoreDiagonal >= scoreLeft && scoreDiagonal >= scoreRight) printf("*\t");
            //else if(scoreRight > scoreLeft) printf("{\t"); else printf("}\t");
            //getchar();
            #endif
            j_prime++;
        }
        #ifdef VERBOSE
        printf("\n");
        getchar();
        #endif
    }
        
    return bc;
}



void backtrackingNW(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, struct cell ** table, char * rec_X, char * rec_Y, struct best_cell * bc, uint64_t * ret_head_x, uint64_t * ret_head_y, int64_t * cell_path_y, uint64_t window_size, BasicAlignment * ba){
    uint64_t curr_x, curr_y, prev_x, prev_y, head_x, head_y;
    int64_t k, j_prime, delta_diff = 0;
    head_x = 2*MAX(Xend, Yend);
    head_y = 2*MAX(Xend, Yend);
    curr_x = bc->c.xpos;
    curr_y = bc->c.ypos;
    #ifdef VERBOSE
    printf("Optimum : %"PRIu64", %"PRIu64"\n", curr_x, curr_y);
    #endif
    prev_x = curr_x;
    prev_y = curr_y;
   
    for(k=Xend-1; k>curr_x; k--) rec_X[head_x--] = '-';
    for(k=Yend-1; k>curr_y; k--) rec_Y[head_y--] = '-';

    j_prime = bc->j_prime;
    unsigned char first_track = 1;
    
    while(curr_x > 0 && curr_y > 0){

        
        if(first_track == 0){
            delta_diff = MAX(1, cell_path_y[prev_x] - (int64_t) window_size) - MAX(1, cell_path_y[curr_x] - (int64_t)window_size); //j-1
            j_prime = j_prime - (int64_t)(prev_y - curr_y) + (int64_t) delta_diff;

            prev_x = curr_x;
            prev_y = curr_y;

            #ifdef VERBOSE
            //printf("Jprime: %"PRId64" :DELTADIF:%"PRId64"\n", j_prime, delta_diff);
            printf("[%c %c]", X[prev_x], Y[prev_y]);
            printf("(%"PRIu64", %"PRIu64") ::: \n", curr_x, curr_y);
            //printf("(%"PRIu64", %"PRIu64") ::: \n", prev_x, prev_y);
            //printf("cellp Prev: %"PRId64" Post: %"PRId64"\n", cell_path_y[prev_x], cell_path_y[curr_x]);
            //printf("the difs? %"PRId64" the other: %"PRId64"\n", MAX(0, cell_path_y[prev_x] - (int64_t) window_size), MAX(0, cell_path_y[curr_x] - (int64_t)window_size));
            getchar();
            #endif

        }

        curr_x = table[prev_x][j_prime].xfrom;
        curr_y = table[prev_x][j_prime].yfrom;
        first_track = 0;
        

        

        if((curr_x == (prev_x - 1)) && (curr_y == (prev_y -1))){
            //Diagonal case
            //printf("DIAG\n");
            rec_X[head_x--] = (char) X[prev_x];
            rec_Y[head_y--] = (char) Y[prev_y];
            ba->length++;
            
        }else if((prev_x - curr_x) > (prev_y - curr_y)){
            //Gap in X
            //printf("Gap X\n");
            for(k=prev_x-1;k>curr_x;k--){
                rec_Y[head_y--] = '-';
                rec_X[head_x--] = (char) X[k];
                ba->length++;
                ba->egaps++;
            }
            ba->igaps += 1;
            ba->egaps--;
        }else{
            //Gap in Y
            //printf("GAP Y\n");
            //10, 0, 401, 281
            for(k=prev_y-1;k>curr_y;k--){
                rec_X[head_x--] = '-';
                rec_Y[head_y--] = (char) Y[k];
                ba->length++;
                ba->egaps++;
            }
            ba->igaps += 1;
            ba->egaps--;
        }
        
    }
    if(curr_x == 0 && curr_y == 0){
        rec_X[head_x--] = (char) X[prev_x];
        rec_Y[head_y--] = (char) Y[prev_y];
        ba->length++;
    }
    
    //printf("curr: %"PRIu64", %"PRIu64"\n", curr_x, curr_y);
    //printf("Heads: %"PRIu64", %"PRIu64"\n", head_x, head_y);
    uint64_t huecos_x = 0, huecos_y = 0;
    for(k=(int64_t)curr_x-1; k>=0; k--){ rec_X[head_x--] = '-'; huecos_x++;}
    for(k=(int64_t)curr_y-1; k>=0; k--){ rec_Y[head_y--] = '-'; huecos_y++;}
    
    if(huecos_x >= huecos_y){
        while(huecos_x > 0) {rec_Y[head_y--] = ' '; huecos_x--;}
    }else{
        while(huecos_y > 0) {rec_X[head_x--] = ' '; huecos_y--;}
    }

    *ret_head_x = head_x;
    *ret_head_y = head_y;
    #ifdef VERBOSE
    printf("hx hy: %"PRIu64", %"PRIu64"\n", head_x, head_y);
    #endif
}
