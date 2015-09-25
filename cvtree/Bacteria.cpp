//
//  Bacteria.cpp
//  cvtree
//
//  Created by Luiz Felix on 20/09/2015.
//  Copyright (c) 2015 Luiz Felix. All rights reserved.
//

#include "Bacteria.h"

const short Bacteria::code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};

const long Bacteria::M2 = (long)(pow(AA_NUMBER, LEN - 2));
const long Bacteria::M1 = M2 * AA_NUMBER;
const long Bacteria::M = M1 * AA_NUMBER;

/* Initialises the arrays on 0, used on computation for each bacteria */
void Bacteria::InitVectors() {
    Bacteria:vector = new long [M];
    second = new long [M1];
    memset(second, 0, M1 * sizeof(long));
    memset(vector, 0, M * sizeof(long));
    memset(one_l, 0, AA_NUMBER * sizeof(long));
    total = 0;
    total_l = 0;
    complement = 0;
}

/* Appends the initial aminoacid sequence to the arrays */
void Bacteria::init_buffer(char* buffer)
{
    complement++;
    indexs = 0;
    for (int i=0; i<LEN-1; i++)
    {
        short enc = encode(buffer[i]);
        one_l[enc]++;
        total_l++;
        indexs = indexs * AA_NUMBER + enc;
    }
    second[indexs]++;
}

/* Append a new aminoacid to the arrays and update signature */
void Bacteria::cont_buffer(char ch) {
    short enc = encode(ch);
    one_l[enc]++;
    total_l++;
    long index = indexs * AA_NUMBER + enc;
    vector[index]++;
    total++;
    indexs = (indexs % M2) * AA_NUMBER + enc;
    second[indexs]++;
}

/* Creates a new bacteria and fills the aminoacids arrays */
Bacteria::Bacteria(char *filename) {
    FILE * bacteria_file = fopen(filename,"r");
    InitVectors();
    
    /* This part is the same */
    char ch;
    while ((ch = fgetc(bacteria_file)) != EOF)
    {
        if (ch == '>')
        {
            while (fgetc(bacteria_file) != '\n'); // skip rest of line
            
            char buffer[LEN-1];
            fread(buffer, sizeof(char), LEN-1, bacteria_file);
            init_buffer(buffer);
        }
        else if (ch != '\n' && ch != '\r')
            cont_buffer(ch);
    }
    
    fclose (bacteria_file);
}

void* threaded_data(void *args) {
    Bacteria::thread_parameters *params = (Bacteria::thread_parameters*)args;
    
    long lower_bound = params->thread_id * params->block_size * Bacteria::M1;
    long upper_bound = (params->thread_id == params->pool_size - 1) ?  Bacteria::M : (params->thread_id + 1) * params->block_size * Bacteria::M1;
    
    long internal_count = 0;
    
    double total_div_2 = params->half_total;
    long i_mod_aa_number = 0;
    long i_div_aa_number = params->thread_id * params->block_size * Bacteria::M2;
    long i_mod_M1 = 0;
    long i_div_M1 = params->thread_id * params->block_size;
    
    for (long i = lower_bound; i < upper_bound; i++) {
        double p1 = params->second_div_total[i_div_aa_number];
        double p2 = params->one_l_div_total[i_mod_aa_number];
        double p3 = params->second_div_total[i_mod_M1];
        double p4 = params->one_l_div_total[i_div_M1];
        double stochastic =  (p1 * p2 + p3 * p4) * total_div_2;
        
        if (i_mod_M1 == Bacteria::M1-1)
        {
            i_mod_M1 = 0;
            i_div_M1++;
        }
        else
            i_mod_M1++;
        
        if (i_mod_aa_number == AA_NUMBER-1)
        {
            i_mod_aa_number = 0;
            i_div_aa_number++;
        }
        else
            i_mod_aa_number++;
        
        if (stochastic > EPSILON)
        {
            params->t[i] = (params->vector[i] - stochastic) / stochastic;
            internal_count++;
        }
        else
            params->t[i] = 0;
    }
    
    // return the amount of values > EPSILON
    pthread_exit((void*)internal_count);
}

/* Does the stochastic >heavy< computation */
void Bacteria::stochastic(pthread_t thread_pool[], int pool_size) {
    
    long total_plus_complement = total + complement;
    
    double one_l_div_total[AA_NUMBER];
    for (int i=0; i<AA_NUMBER; i++)
        one_l_div_total[i] = (double)one_l[i] / total_l;
    
    double* second_div_total = new double[M1];
    for (int i=0; i<M1; i++)
        second_div_total[i] = (double)second[i] / total_plus_complement;
    
    double* t = new double[M];
    
    // dividing the stochastic task
    if (0 > pool_size || pool_size > AA_NUMBER) {
        printf("Amount of threads for bacteria initialisation must be between [0, %d]\n", AA_NUMBER);
        exit(-1);
    }
    
    // enforces joinable threads
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    long big_block = M/M1;  // (k == AA)
    
    thread_parameters params[pool_size];
    
    for (int i = 0; i < pool_size; i++) {
        params[i].thread_id = i;
        params[i].pool_size = pool_size;
        params[i].block_size = big_block / pool_size;
        params[i].half_total = total * 0.5;
        
        // setting up pointers
        params[i].one_l_div_total = one_l_div_total;
        params[i].second_div_total = second_div_total;
        
        params[i].t = t;
        params[i].vector = vector;
        
        int rc = pthread_create(&thread_pool[i], &attr, threaded_data, (void*)&params[i]);
        
        if (rc) {
            printf("Error while creating bacteria thread.\n");
            exit(1);
        }
    }
    
    count = 0;
    for (long i = 0; i < pool_size; i++) {
        void *returnValue;
        pthread_join(thread_pool[i], &returnValue);
        
        count += (long)returnValue;
    }

    delete second_div_total;
    delete vector;
    delete second;
    
    tv = new double[count];
    ti = new long[count];
    
    int pos = 0;
    double sum = 0;
    for (long i=0; i<M; i++)
    {
        if (t[i] != 0)
        {
            tv[pos] = t[i];
            ti[pos] = i;
            pos++;
            sum += t[i];
        }
    }
    
    delete t;
}
