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

#define encode(ch)		code[ch-'A']
#define LEN				6
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

    static const long M, M1, M2;

    void InitVectors();
    void init_buffer(char* buffer);
    void cont_buffer(char ch);
    
public:
    /* Fields used on the stochastic method */
    long count;
    double* tv;
    long *ti;
    
    Bacteria(char *name);
};


#endif /* defined(__cvtree__Bacteria__) */
