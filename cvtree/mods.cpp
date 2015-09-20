/*
 This program is the same as improved.cpp (supplied), but it times only CPU time.
 The results are stored in a temporary array, allowing to remove the IO time from
 the profiling.
 20/09/2015
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "Bacteria.h"

//#define min(a, b) ((a) < (b)) ? (a) : (b)

namespace modifications {
    
    //same things
    int number_bacterias;
    
    char** bacteria_name;

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
    time_t t1 = time(NULL);
    
    ReadInputFile(argv[1]);
    
    CompareAllBacteria();
    
    time_t t2 = time(NULL);
    
    printf("time elapsed: %ld seconds\n", t2 - t1);
    
    show_similarities();
    
    return 0;
}