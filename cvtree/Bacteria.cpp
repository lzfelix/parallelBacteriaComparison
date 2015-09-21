//
//  Bacteria.cpp
//  cvtree
//
//  Created by Luiz Felix on 20/09/2015.
//  Copyright (c) 2015 Luiz Felix. All rights reserved.
//

#include "Bacteria.h"

const short Bacteria::code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};

//void Init()
//{
//    M2 = 1;
//    
//    for (int i=0; i<LEN-2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
//        M2 *= AA_NUMBER;
//    
//    M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
//    M  = M1 * AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
//}

const long Bacteria::M2 = (long)(pow(AA_NUMBER, LEN - 2));
const long Bacteria::M1 = M2 * AA_NUMBER;
const long Bacteria::M = M1 * AA_NUMBER;

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

void Bacteria::stochastic() {
    // keeping some calculations stored
    long total_plus_complement = total + complement;
    double total_div_2 = total * 0.5;
    int i_mod_aa_number = 0;
    int i_div_aa_number = 0;
    long i_mod_M1 = 0;
    long i_div_M1 = 0;
    
    double one_l_div_total[AA_NUMBER];
    for (int i=0; i<AA_NUMBER; i++)
        one_l_div_total[i] = (double)one_l[i] / total_l;
    
    double* second_div_total = new double[M1];
    for (int i=0; i<M1; i++)
        second_div_total[i] = (double)second[i] / total_plus_complement;
    
    // this is the stochastic processing for comparison. Computes everything and
    // stores. Need to understand better.
    count = 0;
    double* t = new double[M];
    
    for(long i=0; i<M; i++)
    {
        double p1 = second_div_total[i_div_aa_number];
        double p2 = one_l_div_total[i_mod_aa_number];
        double p3 = second_div_total[i_mod_M1];
        double p4 = one_l_div_total[i_div_M1];
        double stochastic =  (p1 * p2 + p3 * p4) * total_div_2;
        
        if (i_mod_aa_number == AA_NUMBER-1)
        {
            i_mod_aa_number = 0;
            i_div_aa_number++;
        }
        else
            i_mod_aa_number++;
        
        if (i_mod_M1 == M1-1)
        {
            i_mod_M1 = 0;
            i_div_M1++;
        }
        else
            i_mod_M1++;
        
        if (stochastic > EPSILON)
        {
            t[i] = (vector[i] - stochastic) / stochastic;
            count++;
        }
        else
            t[i] = 0;
    }
    
    delete second_div_total;
    delete vector;
    delete second;
    
    tv = new double[count];
    ti = new long[count];
    
    int pos = 0;
    for (long i=0; i<M; i++)
    {
        if (t[i] != 0)
        {
            tv[pos] = t[i];
            ti[pos] = i;
            pos++;
        }
    }
    delete t;
}
