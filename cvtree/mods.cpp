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

#define MIN(a, b) ((a) < (b)) ? (a) : (b)

#define NUM_THREADS 1

namespace modifications {
    
    struct parameters {
        int thread_id;
        Bacteria **bacterias;
    };
    
    int number_bacterias;       //amount of bacterias to read
    char** bacteria_name;       //array of bacteria names

    //hardcoded amount of bacteria
    #define BACTERIA_AMOUNT 41
    double similarity_table[BACTERIA_AMOUNT][BACTERIA_AMOUNT];
    
    //path the the folder with the bacteria files
    #define FOLDER_NAME     "data/"
    
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
    
    
    /* Threaded function spawned by [threaded_bacteria_comparison]. Populates the bacteria array */
    void *create_bacterias(void *args) {
        parameters *params = (parameters*)args;
        
        // calculating bounds and workload
        int block_size = (number_bacterias) / NUM_THREADS;
        int lower_bound = params->thread_id * block_size;
        int upper_bound = 0;
        
        if (params->thread_id == NUM_THREADS - 1)
            upper_bound = number_bacterias;
        else
            upper_bound = lower_bound + block_size;
    

        for (int i = lower_bound; i < upper_bound; i++) {
            params->bacterias[i] = new Bacteria(bacteria_name[i]);
        }

        
        printf("Thread %d done.\n", params->thread_id);
        
        return 0;
    }
    
    /* Load the files retrieved from the list parallelly. */
    void threaded_bacteria_comparison() {
        
        // creates the bacterias array
        Bacteria **bacterias = new Bacteria*[number_bacterias];
        
        pthread_t workers[NUM_THREADS];
        parameters params[NUM_THREADS];
        
        // enforces joinable threads
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        
        for (long i = 0; i < NUM_THREADS; i++) {
            
            // the thread parameters
            params[i].thread_id = i;
            params[i].bacterias = bacterias;
            
            int rc = pthread_create(&workers[i], &attr, create_bacterias, (void*)&params[i]);
            
            if (rc)
                terminate_program("Error while creating threads.");
        }
        
        // Joines the thread without caring for the returned value
        for (long i = 0; i < NUM_THREADS; i++)
            pthread_join(workers[i], NULL);
    }
    
    /* Slightly changed */
    void CompareAllBacteria()
    {
        //creates an array of bacterias
        Bacteria** b = new Bacteria*[number_bacterias];
        
        for(int i=0; i<number_bacterias; i++)
        {
            b[i] = new Bacteria(bacteria_name[i]);
        }
        
//        performs the comparison
        for(int i=0; i<number_bacterias-1; i++)
            for(int j=i+1; j<number_bacterias; j++)
                similarity_table[i][j] = CompareBacteria(b[i], b[j]);
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

//    CompareAllBacteria();
    
    gettimeofday(&t_end, NULL);
    
    //calculating time
    elapsed_time = (t_end.tv_sec - t_init.tv_sec) * 1000.0;       //transforming to ms
    elapsed_time += (t_end.tv_usec - t_init.tv_usec) / 1000.0;    //moving to the right place
    
    printf("\nIn %f ms.\n", elapsed_time);
    
//    show_similarities();
    
    return 0;
}