#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <sys/time.h>       // timming

namespace improved {
    //same things
    int number_bacteria;
    char** bacteria_name;
    long M, M1, M2;
    short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
    #define encode(ch)		code[ch-'A']
    #define LEN				6
    #define AA_NUMBER		20
    #define	EPSILON			1e-010
    
    //path the the folder with the bacteria files
    #define FOLDER_NAME     "data/"

    /*
     The same, inits the constants.
    */
    void Init()
    {
        M2 = 1;
        for (int i=0; i<LEN-2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
            M2 *= AA_NUMBER; 
        M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
        M  = M1 *AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
    }

    class Bacteria
    {
    private:
        /*
            Same thing, but now vector is also private.
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

        /* Same thing */
        void InitVectors()
        {
            vector = new long [M];
            second = new long [M1];
            memset(vector, 0, M * sizeof(long));
            memset(second, 0, M1 * sizeof(long));
            memset(one_l, 0, AA_NUMBER * sizeof(long));
            total = 0;
            total_l = 0;
            complement = 0;
        }

        /* Same thing */
        void init_buffer(char* buffer)
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

        /* Same thing */
        void cont_buffer(char ch)
        {
            short enc = encode(ch);
            one_l[enc]++;
            total_l++;
            long index = indexs * AA_NUMBER + enc;
            vector[index]++;
            total++;
            indexs = (indexs % M2) * AA_NUMBER + enc;
            second[indexs]++;
        }

    public:
        /* These guys are new */
        long count;
        double* tv;
        long *ti;

        /* Changed A LOT */
        Bacteria(char* filename)
        {
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

            fclose (bacteria_file);
        }
    };
    //-- class end

    /* Same thing */
    void ReadInputFile(char* input_name)
    {
        FILE* input_file = fopen(input_name,"r");
        fscanf(input_file,"%d",&number_bacteria);
        bacteria_name = new char*[number_bacteria];

        for(long i=0;i<number_bacteria;i++)
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
            long n1 = b1->ti[p1];
            double t1 = b1->tv[p1++];
            vector_len1 += (t1 * t1);
        }
        
        while (p2 < b2->count)
        {
            long n2 = b2->ti[p2];
            double t2 = b2->tv[p2++];
            vector_len2 += (t2 * t2);
        }
        
        //face here, duh

        return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
    }

    /* Slightly changed */
    void CompareAllBacteria()
    {
        //This can be interleaved.
        
        //loads all bacterias in an array
        Bacteria** b = new Bacteria*[number_bacteria];
        for(int i=0; i<number_bacteria; i++)
        {
            b[i] = new Bacteria(bacteria_name[i]);
        }

        //performs the comparison
        for(int i=0; i<number_bacteria-1; i++)
            for(int j=i+1; j<number_bacteria; j++)
            {
                printf("%2d %2d -> ", i, j);
                double correlation = CompareBacteria(b[i], b[j]);
                printf("%.20lf\n", correlation);
            }
    }
}

//using namespace improved;
//
//int main(int argc,char * argv[])
//{
//    struct timeval t_init, t_end;
//    double elapsed_time;
//    
//    gettimeofday(&t_init, NULL);
//
//	Init();
//	ReadInputFile(argv[1]);
//	CompareAllBacteria();
//    
//    gettimeofday(&t_end, NULL);
//
//    //calculating time
//    elapsed_time = (t_end.tv_sec - t_init.tv_sec) * 1000.0;       //transforming to ms
//    elapsed_time += (t_end.tv_usec - t_init.tv_usec) / 1000.0;    //moving to the right place
//    
//    printf("In %f ms.\n", elapsed_time);
//    
//    return 0;
}