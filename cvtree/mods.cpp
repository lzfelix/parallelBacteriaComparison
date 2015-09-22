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

#define NUM_THREADS 4

//path the the folder with the bacteria files
#define FOLDER_NAME "data/"

namespace modifications {
    
    struct parameters {
        int thread_id;
        int bacterias_read;                // the amount of already bacteria files read
        Bacteria **bacterias;             // pointer to the array of bacterias
        pthread_mutex_t mutex;      // the mutex to protect [bacterias_read]
        pthread_cond_t  signal;        // used to wake the thread that runs ::stochastic() over each bacteria on the thread block
        
        // bounding variables
        int lower_bound;
        int upper_bound;
        int block_size;
    };
    
    int number_bacterias;       //amount of bacterias to read
    char** bacteria_name;       //array of bacteria names

    //hardcoded amount of bacteria
    #define BACTERIA_AMOUNT 41
    double similarity_table[BACTERIA_AMOUNT][BACTERIA_AMOUNT];
    
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
        
        //These two loops can be paralellised
        
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
        
        //face here, duh
        
        return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
    }
    
    /* Displays the error message [msg] and exists with status code 1 */
    void terminate_program(const char* msg) {
        printf("%s\n", msg);
        exit(1);
    }
    
    /* For each block, finds the lower, upper bounds and block size. Call this just after having parsed the list file */
    /* inline */ void calculate_block_parameters(parameters *p) {
        int id = p->thread_id;
        
        p->block_size = (number_bacterias) / NUM_THREADS;
        p->lower_bound = p->thread_id * p->block_size;
        
        if (p->thread_id == NUM_THREADS - 1)
            p->upper_bound = number_bacterias;
        else
            p->upper_bound = p->lower_bound + p->block_size;
        
    }
    
    /* Threaded function spawned by [threaded_bacteria_comparison]. Populates the bacteria array */
    void *create_bacterias(void *args) {
        parameters *params = (parameters*)args;
    
        //no bacterias were read so far
        params->bacterias_read = 0;
        
        int deltaSize = params->block_size * 0.25;
        int max = params->upper_bound - params->lower_bound;

        for (int i = params->lower_bound; i < params->upper_bound; i++) {
            params->bacterias[i] = new Bacteria(bacteria_name[i]);
            params->bacterias[i]->stochastic();
            
            pthread_mutex_lock(&(params->mutex));
            params->bacterias_read++;

            if (params->bacterias_read > 0 && (params->bacterias_read % deltaSize == 0 || params->bacterias_read == max) ){
                pthread_cond_signal(&(params->signal));
            }
            
            pthread_mutex_unlock(&(params->mutex));
        }
        
        return 0;
    }
    

    void *analyse_bacteria(void *args) {
        parameters *params = (parameters*)args;
    
        int block_index = params->lower_bound;
        int upper_bound = params->upper_bound;
        
        pthread_mutex_lock(&(params->mutex));
        
        while (block_index < upper_bound) {
            // complying with pthreads spec.
            while (params->bacterias[block_index] == NULL) {
                pthread_cond_wait(&(params->signal), &(params->mutex));
//                printf("Locked %d, %d\n", relative_index, params->bacterias_read);
            }
            
//            printf("#%d Analyser analysing %d\n", params->thread_id, block_index);
            params->bacterias[block_index]->stochastic();
            block_index++;
        }
        
        pthread_mutex_unlock(&(params->mutex));
        
        
        return 0;
    }
    
    /* Load the files retrieved from the list parallelly. */
    void threaded_bacteria_comparison() {
        
        // creates the bacterias array
        Bacteria **bacterias = new Bacteria*[number_bacterias];
        
        pthread_t workers[NUM_THREADS];
        pthread_t analysers[NUM_THREADS];
        parameters params[NUM_THREADS];
        
        // enforces joinable threads
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        
        // spawning thread responsible for creating bacterias
        for (long i = 0; i < NUM_THREADS; i++) {
            // the thread parameters
            params[i].thread_id = i;
            params[i].bacterias = bacterias;
            
            // these parameters are only used analyse_bacteria()
            pthread_mutex_init(&params[i].mutex, NULL);
            pthread_cond_init(&params[i].signal, NULL);
            
            calculate_block_parameters(&params[i]);
            
            int rc = pthread_create(&workers[i], &attr, create_bacterias, (void*)&params[i]);
            
            if (rc)
                terminate_program("Error while creating threads.");
        }
        
        // spawns thread that analyse each bacteria
        for (int i = 0; i < NUM_THREADS; i++) {
            pthread_create(&analysers[i], &attr, analyse_bacteria, (void*)&params[i]);
        }
        
        // joins the worker threads without caring for the returned value
        for (long i = 0; i < NUM_THREADS; i++)
            pthread_join(workers[i], NULL);
        
        // joins the analysers
        for (long i = 0; i < NUM_THREADS; i++)
            pthread_join(analysers[i], NULL);
        
////        perform comaparison
        for(int i=0; i<number_bacterias-1; i++)
            for(int j=i+1; j<number_bacterias; j++)
                similarity_table[i][j] = CompareBacteria(bacterias[i], bacterias[j]);
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
    threaded_bacteria_comparison();
    
    gettimeofday(&t_end, NULL);
    
    //calculating time
    elapsed_time = (t_end.tv_sec - t_init.tv_sec) * 1000.0;       //transforming to ms
    elapsed_time += (t_end.tv_usec - t_init.tv_usec) / 1000.0;    //moving to the right place
    
    printf("\nIn %f ms.\n", elapsed_time);
    
//    show_similarities();
    
    return 0;
}