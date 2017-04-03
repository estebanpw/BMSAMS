

/*
    Piles up a series of chars so that we can produce scores
*/

void pile_chars_up(char ** list_seq_X, char ** list_seq_Y, char * list_piled_X, char * list_piled_Y, uint64_t nx, uint64_t ny, uint64_t posx, uint64_t posy);
/*
    Compares nucleotides of a column in a multiple sequence alignment NW
*/
int64_t compare_piled_up_chars(char * list_piled_X, char * list_piled_Y, uint64_t nx, uint64_t ny);

/*
    Nucleotides matching function
*/
inline int64_t compare_letters(char a, char b);

/*
    Gets the maximum length of a group of sequences
*/
uint64_t get_max_length_of_sequences(char ** a, uint64_t n_x);
/*
    Builds a cellpath of diagonal steps for a bounded NW when no frags are available
*/
void build_unanchored_alignment(int64_t * cell_path_y, char ** seq_piles_up_X, char ** seq_piles_up_Y, uint64_t n_x, uint64_t n_y);

/*
    Performs NW and backtracking to recover alignment
*/
All_seqs build_multiple_alignment(char ** reconstruct_X, char ** reconstruct_Y, char ** my_x, char ** my_y, uint64_t nx, uint64_t ny, struct cell ** table, struct positioned_cell * mc, char * writing_buffer_alignment, uint64_t xlen, uint64_t ylen, int64_t * cell_path_y, long double * window, int64_t iGap, int64_t eGap, BasicAlignment * ba, char * aux);


/*
    Computes the cell path for the y points given incremental x
    Only add +- window size to each to know which path to go through
*/
void calculate_y_cell_path(Point p0, Point p1, Point p2, Point p3, int64_t * cell_path_y);
/*
    Calculates NW table with two rows and stores a cellpath of scores, identities, gaps and starting and ending positions
*/
struct best_cell NW(char ** X, uint64_t Xstart, uint64_t Xend, char ** Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell ** table, struct positioned_cell * mc, int show, int64_t * cell_path_y, long double * window, uint64_t * current_window_size, uint64_t nx, uint64_t ny);

/*
    Computes the alignment given a NW table
*/
void backtrackingNW(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, struct cell ** table, char * rec_X, char * rec_Y, struct best_cell * bc, uint64_t * ret_head_x, uint64_t * ret_head_y, int64_t * cell_path_y, uint64_t window_size, BasicAlignment * ba);

