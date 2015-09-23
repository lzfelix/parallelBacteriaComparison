/*
 This program is the same as improved.cpp (supplied), but it times only CPU time.
 The results are stored in a temporary array, allowing to remove the IO time from
 the profiling.
 20/09/2015
 */

#include <cstdio>
#include <cstdlib>

#include "Bacteria.h"
#include <pthread.h>
#include <sys/time.h>       // timming

#define NUM_THREADS 5

//path the the folder with the bacteria files
#define FOLDER_NAME "data/"

namespace modifications {
    
    struct parameters {
        int thread_id;
        int bacterias_read;                // the amount of already bacteria files read
        Bacteria **bacterias;             // pointer to the array of bacterias
        
        int lower_bound;
        int upper_bound;
    };

    int number_bacterias;       //amount of bacterias to read
    char** bacteria_name;       //array of bacteria names
    
    //these tables are used to distribute the comparison workload
    //(which is a upper triangular matrix with main diagonal = 0)
    //evenly among threads
    int *triangle_table;
    int *triangle_bounds;

    //hardcoded amount of bacteria
    #define BACTERIA_AMOUNT 41
    double similarity_table[BACTERIA_AMOUNT][BACTERIA_AMOUNT];
    
    /* Computes the first and last value of each line on the triangular matrix */
    void create_triangle_table() {
        triangle_table = new int[number_bacterias - 1];
        triangle_bounds = new int[number_bacterias];
        
        int seed = number_bacterias - 2;
        int seed_sum = seed;
        int lower_bound = 0;
        
        for (int i = 0; i < number_bacterias - 1; i++) {
            triangle_table[i] = seed;
            triangle_bounds[i] = lower_bound;
            
            lower_bound = seed + 1;
            seed += seed_sum--;
        }
        
        triangle_bounds[number_bacterias - 1] = 10000;
    }
    
    /* finds the next element to search on the comparison table */
    inline int pseudo_binary_search(int key, int a, int b) {
        int middle = (a + b) / 2;
        
        if (key > triangle_table[number_bacterias - 2] || key < 0)
            return -1;
        
        if (key <= triangle_table[middle]) {
            if (middle - 1 < 0) return 0;
            
            if (key > triangle_table[middle - 1])
                return middle;
            else
                return pseudo_binary_search(key, a, middle);
        }
        else
            return pseudo_binary_search(key, middle, b);
    }
    
    /* Same thing */
    void ReadInputFile(char* input_name)
    {
        FILE* input_file = fopen(input_name,"r");
        fscanf(input_file,"%d",&number_bacterias);
        bacteria_name = new char*[number_bacterias];
        
        for(long i=0;i<number_bacterias;i++)
        {
            //contains the full path
            bacteria_name[i] = new char[26];
            
            //the real bacteria name
            char bacteriaName[20];
            
            fscanf(input_file, "%s", bacteriaName);
            
            //-> path/bacteria_name.faa
            strcpy(bacteria_name[i], FOLDER_NAME);
            strcat(bacteria_name[i], bacteriaName);
            strcat(bacteria_name[i],".faa");
        }
        fclose(input_file);
    }
    
    /* Different */
    double CompareBacteria(Bacteria* b1, Bacteria* b2)
    {
        double correlation = 0;
        double vector_len1=0;
        double vector_len2=0;
        long p1 = 0;
        long p2 = 0;
        
        while (p1 < b1->count && p2 < b2->count)
        {
            long n1 = b1->ti[p1];
            long n2 = b2->ti[p2];
            if (n1 < n2)
            {
                double t1 = b1->tv[p1];
                vector_len1 += (t1 * t1);
                p1++;
            }
            else if (n2 < n1)
            {
                double t2 = b2->tv[p2];
                p2++;
                vector_len2 += (t2 * t2);
            }
            else
            {
                double t1 = b1->tv[p1++];
                double t2 = b2->tv[p2++];
                vector_len1 += (t1 * t1);
                vector_len2 += (t2 * t2);
                correlation += t1 * t2;
            }
        }
        
        while (p1 < b1->count)
        {
            double t1 = b1->tv[p1++];
            vector_len1 += (t1 * t1);
        }
        
        while (p2 < b2->count)
        {
            double t2 = b2->tv[p2++];
            vector_len2 += (t2 * t2);
        }
        
        return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
    }
    
    /* Displays the error message [msg] and exists with status code 1 */
    void terminate_program(const char* msg) {
        printf("%s\n", msg);
        exit(1);
    }
    
    /* Threaded function spawned by [threaded_bacteria_comparison]. Populates the bacteria array */
    void *create_bacterias(void *args) {
        parameters *params = (parameters*)args;
        
        for (int i = params->lower_bound; i < params->upper_bound; i++) {
            params->bacterias[i] = new Bacteria(bacteria_name[i]);
            params->bacterias[i]->stochastic();
        }
        
        return 0;
    }
    
    /* Compare a single bacteria with all others evenly */
    void *threaded_compare_bacterias(void *args) {
        parameters *params = (parameters*)args;
        
        int r_x = 0;
        int r_y = 0;
        int last_y = 0;
        int max_x = -1;
     
         for (int i = params->lower_bound; i < params->upper_bound; i++) {
             if ( i > max_x) {
                 r_y = pseudo_binary_search(i, 0, number_bacterias - 1);
                 last_y = r_y;
                 max_x = triangle_bounds[r_y + 1] - 1;
             }
             else
                 r_y = last_y;
             
             r_x = r_y + 1 + i - triangle_bounds[r_y];
         
             similarity_table[r_y][r_x] = CompareBacteria(params->bacterias[r_y], params->bacterias[r_x]);
         }
    
        return 0;
    }
    
    /* Load the files retrieved from the list parallelly. */
    void threaded_bacteria_creation() {
        
        // creates the bacterias array
        Bacteria **bacterias = new Bacteria*[number_bacterias];
        
        pthread_t workers[NUM_THREADS];
        parameters params[NUM_THREADS];
        
        // enforces joinable threads
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        
        int block_size = number_bacterias / NUM_THREADS;
        
        // spawns threads responsible for creating bacterias
        for (int i = 0; i < NUM_THREADS; i++) {
            
            // the thread parameters
            params[i].thread_id = i;
            params[i].bacterias = bacterias;
            
            // setting the bounds
            params[i].lower_bound = i * block_size;
            
            if (params[i].thread_id == NUM_THREADS - 1)
                params[i].upper_bound = number_bacterias;
            else
                params[i].upper_bound = params[i].lower_bound + block_size;
            
            // creating the thread
            int rc = pthread_create(&workers[i], &attr, create_bacterias, (void*)&params[i]);
            
            if (rc)
                terminate_program("Error while creating threads.");
        }
        
        // joins the worker threads without caring for the returned value
        for (long i = 0; i < NUM_THREADS; i++)
            pthread_join(workers[i], NULL);
        
        // create threads to perform the comparison
        int max_comparisons = number_bacterias * 0.5 * (number_bacterias - 1);
        block_size = max_comparisons / NUM_THREADS;
        
        for (int i = 0; i < NUM_THREADS; i++) {
            params[i].lower_bound = block_size * i;
            
            if (params[i].thread_id == NUM_THREADS - 1)
                params[i].upper_bound = max_comparisons;
            else
                params[i].upper_bound = params[i].lower_bound + block_size;
            
            int rc = pthread_create(&workers[i], &attr, threaded_compare_bacterias, (void*)&params[i]);
            
            if (rc)
                terminate_program("Error while creating comparison threads.\n");
        }
        
        // join the comparison threads
        for (long i = 0; i < NUM_THREADS; i++)
            pthread_join(workers[i], NULL);
    }
    
    void show_similarities() {
        for (int i = 0; i < BACTERIA_AMOUNT; i++)
            for (int j = i + 1; j < BACTERIA_AMOUNT ; j++)
                printf("%2d %2d -> %.20lf\n", i, j, similarity_table[i][j]);
    }
}

using namespace modifications;

int main(int argc,char * argv[])
{
    struct timeval t_init, t_end;
    double elapsed_time;
    
    gettimeofday(&t_init, NULL);

    ReadInputFile(argv[1]);
    
    create_triangle_table();
    threaded_bacteria_creation();
    
    gettimeofday(&t_end, NULL);
    
    //calculating time
    elapsed_time = (t_end.tv_sec - t_init.tv_sec) * 1000.0;       //transforming to ms
    elapsed_time += (t_end.tv_usec - t_init.tv_usec) / 1000.0;    //moving to the right place
    
    printf("In %f ms.\n", elapsed_time);
    
    show_similarities();
    
    return 0;
}