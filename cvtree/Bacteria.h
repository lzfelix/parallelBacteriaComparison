//
//  Bacteria.h
//  An header representation of the Bacteria class.
//
//  Created by Luiz Felix on 20/09/2015.
//  Copyright (c) 2015 Luiz Felix. All rights reserved.
//

#ifndef __cvtree__Bacteria__
#define __cvtree__Bacteria__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <pthread.h>

#define encode(ch)		code[ch-'A']
#define LEN                  6
#define AA_NUMBER		20
#define	EPSILON			1e-010

//making methods public

class Bacteria {
private:
    /*
     Information about the fields
     - vector
     - second     : histogram of signatures
     - one_l      : histogram of aminoacids frequency
     - indexs     : the current signature
     - total      : total amount of aminoacids, not in the initial signature
     - total_l    : the real total of aminoacids
     - complement : the amount of signatures (?)
     */
    
    long* vector;
    long* second;
    long one_l[AA_NUMBER];
    long indexs;
    long total;
    long total_l;
    long complement;
    
    static const short code[27];

    void InitVectors();
    void init_buffer(char* buffer);
    void cont_buffer(char ch);
    
    // thread function
//    void* threaded_data(void *args);
    
public:
    struct thread_parameters {
        int thread_id;
        int pool_size;
        long block_size;
        double half_total;
        
        //pointers to arrays
        double *one_l_div_total;
        double *second_div_total;
        double *t;
        long *vector;
    };
    
        static const long M, M1, M2;
    
    /* Fields used on the stochastic method */
    long count;
    double* tv;
    long *ti;
    
    friend void *threaded_data(void *args);
    Bacteria(char *name);
    
    //previously this method was inside the constructor
    void stochastic(pthread_t thread_pool[], int size);
};


#endif /* defined(__cvtree__Bacteria__) */
