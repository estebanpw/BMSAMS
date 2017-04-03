#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include "structs.h"

void terror(char *s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

Queue * generate_queue(Head * queue_head, uint64_t t_reads, uint64_t n_threads, uint64_t levels){
    uint64_t i, j;
    uint64_t reads_per_thread;
    if(levels > t_reads) levels = t_reads;
    uint64_t pieces = t_reads/levels;
    uint64_t from, to, t_queues = 0, current_queue = 0;
    for(i=0;i<levels;i++) t_queues += ((i+1)*n_threads);
    Queue * queues = (Queue *) malloc(t_queues * sizeof(Queue));
    if(queues == NULL) terror("Could not allocate queue tasks");
    queue_head->head = &queues[0];
    
    for(i=0;i<levels;i++){

        //reads_per_thread = (uint64_t) (floorl((long double) pieces / (long double) ((i+1)*n_threads)));
        reads_per_thread = (uint64_t) (ceill((long double) pieces / (long double) ((i+1)*n_threads)));
        

        for(j=0;j<(i+1)*n_threads;j++){
            from = j * reads_per_thread + (pieces*i);
            to = (j + 1) * reads_per_thread + (pieces*i);
            
            if(j == (i+1)*n_threads - 1) to = pieces*(i+1);


            if(i == levels - 1 && j == (i+1)*n_threads - 1){
                //If its the last 
                queues[current_queue].next = NULL;
            }else{
                //Else add it to the queue
                queues[current_queue].next = &queues[current_queue+1];
            }
            

            queues[current_queue].r1 = from;
            queues[current_queue].r2 = to;
            current_queue++;
            //printf("current_piece: %"PRIu64"-%"PRIu64" diff: %"PRIu64"\n", from, to, to - from);

        }

    }
    //printf("TREADS was %"PRIu64"\n", t_reads);    
    return &queues[0];
}

void print_queue(Queue * q){
    fprintf(stdout, "Task: %"PRIu64"-%"PRIu64"\n", q->r1, q->r2);
}

Queue * get_task_from_queue(Head * queue_head, pthread_mutex_t * lock){
    pthread_mutex_lock(lock);

    Queue * ptr = queue_head->head;
    if(queue_head->head != NULL) queue_head->head = queue_head->head->next;
    if(ptr != NULL){ printf("Taking "); print_queue(ptr); }

    pthread_mutex_unlock(lock);


    return ptr;
}
