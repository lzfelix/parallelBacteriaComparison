#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <time.h>
#include <math.h>

namespace tidy {
    //amount of bacterias in a file
    int number_bacteria;

    //name of the current bacteria
    char** bacteria_name;
    long M, M1, M2;

    // the aminoacids codes
    short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};

    //path the the folder with the bacteria files
    #define FOLDER_NAME     "data/"

    //transform number stored in a char as a int
    #define encode(ch)		code[ch-'A']

    //The size of the buffer, used when reading the bacteria file
    #define LEN				6

    //the amount of aminoacids
    #define AA_NUMBER		20

    #define	EPSILON			1e-010

    /**
     Initialises constants M, M1 and M2
    */
    void Init()
    {
        M2 = 1;
        
        for (int i=0; i<LEN-2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
            M2 *= AA_NUMBER;
        
        M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
        M  = M1 * AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
    }


    class Bacteria
    {
        private:
            // the histogram of signatures
            long* second;
        
            // the histogram of aminoacids frequency
            long one_l[AA_NUMBER];
            long indexs;
        
            // who are these guys?
            long total;
        
            // total amount of aminoacids
            long total_l;
            long complement;
        
            /*
             Dynamically creates vector, second, fill them with 0 and set the private long variables to 0.
            */
            void InitVectors()
            {
                // don't know who these guys are (M, M1 length) -> ??????
                vector = new long [M];  //this is public
                second = new long [M1]; //this is private
                
                // fill them with 0s
                memset(vector, 0, M * sizeof(long));
                memset(second, 0, M1 * sizeof(long));
                memset(one_l, 0, AA_NUMBER * sizeof(long));
                
                // ????????????
                total = 0;
                total_l = 0;
                complement = 0;
            }
        
            /*
             For each 6 aminoacids, updates the histogram, total amount of aminoacids and total count.
            */
            void init_buffer(char* buffer)
            {
                complement++;
                indexs = 0;
                for (int i=0; i<LEN-1; i++)
                {
                    //find the number for each aminoacid letter
                    short enc = encode(buffer[i]);
                    
                    //is an histogram of the aminoacids on the bacteria
                    one_l[enc]++;
                    
                    //total amount of aminoacids (?)
                    total_l++;
                    
                    //signature = this encode * 20 + the previous ones. (Stochastic feels).
                    //the signature of the first aminoacid is kept
                    indexs = indexs * AA_NUMBER + enc;
                }
                
                //an histogram of signatures
                second[indexs]++;
            }
            
            void cont_buffer(char ch)
            {
                //update the histograms and total amount of aminoacids with this new one
                short enc = encode(ch);
                one_l[enc]++;
                total_l++;
                
                //calculates the local signature based on the global one.
                long index = indexs * AA_NUMBER + enc;
                
                //histogram of local signatures
                vector[index]++;
                
                //total amount of new aminoacids? (slides?)
                total++;
                
                //updating the global signature
                indexs = (indexs % M2) * AA_NUMBER + enc;
                
                //updating the histogram of global signatures
                second[indexs]++;
            }
            
    public:
            //histogram of local signatures
            long* vector;
            
            /*
             The constructor. Takes a bacteria file name and read such file.
            */
            Bacteria(char* filename)
            {
                FILE * bacteria_file = fopen(filename,"r");
                
                //create the private fields and set them to 0.
                InitVectors();
                
                /*
                 Parsing the file =D
                */
                
                char ch;
                
                //reads each entry on the bacteria file (delimited by >)
                while ((ch = fgetc(bacteria_file)) != EOF)
                {
                    if (ch == '>')
                    {
                        //skips the first line, as it contains meta docs
                        while (fgetc(bacteria_file) != '\n'); // skip the first line, which contains reference data
                        
                        //reads the first 6 aminoacids to find the global signature
                        char buffer[LEN-1];
                        fread(buffer, sizeof(char), LEN-1, bacteria_file);
                        init_buffer(buffer);
                    }
                    
                    //for all others update signature and does other stuff
                    else if (ch != '\n' && ch != '\r')
                        cont_buffer(ch);
                }
                
                fclose (bacteria_file);
            }
            
            /*
             Deletes dynamically allocated vectors that holds both local and global signature.
            */
            ~Bacteria()
            {
                delete vector;
                delete second;
            }
            
            /*
             Can this be stored in memory?
            */
            double stochastic_compute(long i)
            {
                double p1 = (double) second[i / AA_NUMBER] / (total + complement);
                double p2 = (double) one_l[i % AA_NUMBER] / total_l;
                double p3 = (double) second[i % M1] / (total + complement);
                double p4 = (double) one_l[i / M1] / total_l;
                return total * (p1*p2 + p3*p4) / 2;
            }
    };

    /*
     Finds the resemblance of 2 bacterias
    */
    double CompareBacteria(Bacteria* b1, Bacteria* b2)
    {
        double correlation = 0;
        double vector_len1=0;
        double vector_len2=0;
        
        for(long i=0; i<M; i++)
        {
            double stochastic1 = b1->stochastic_compute(i);
            double t1;
            
            if (stochastic1 > EPSILON)
                t1 = (b1->vector[i] - stochastic1) / stochastic1;
            else
                t1=0;
            
            vector_len1 += (t1 * t1);
            
            double stochastic2 = b2->stochastic_compute(i);
            double t2;
            if (stochastic2 > EPSILON)
                t2 = (b2->vector[i] - stochastic2) / stochastic2;
            else
                t2 = 0;
            vector_len2 += (t2 * t2);
            
            correlation = correlation + t1 * t2;
        }
        
        return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
    }


    /*
     //n^2 bacteria comparison. Maybe this is parallelisable
    */
    void CompareAllBacteria()
    {
        for(int i=0; i<number_bacteria-1; i++)
        {
            //create a new bacteria. The comparisor
            Bacteria* b1 = new Bacteria(bacteria_name[i]);
            
            for(int j=i+1; j<number_bacteria; j++)
            {
                //create another bacteria. The compared
                Bacteria* b2 = new Bacteria(bacteria_name[j]);
                
                //compare the bacterias. This is the [[source of heat]]
                double correlation = CompareBacteria(b1, b2);
                printf("%03d %03d -> %.10lf\n", i, j, correlation);
                
                //delete the useless bacteria
                delete b2;
            }
            delete b1;
        }
    }


    /*
     Read bacteria name in files. Changed this method
    */
    void ReadInputFile(char* input_name)
    {
        FILE* input_file = fopen(input_name,"r");
        
        //first element is the amount of bacterias in the file
        fscanf(input_file,"%d",&number_bacteria);
        
        //array of strings containning the name of all bacterias
        bacteria_name = new char*[number_bacteria];
        
        //read each bacteria name and append the .faa extension -- I changed the code here
        for(long i=0;i<number_bacteria;i++)
        {
            //I changed from 20 to 26 -- this is the name with path
            bacteria_name[i] = new char[26];
            
            //add the folder that contains the bacteria data
            strcpy(bacteria_name[i], FOLDER_NAME);
            
            //the bacteria in file name
            char bacteriaName[20];
            
            //read name from file
            fscanf(input_file, "%s", bacteriaName);
            
            //add all togheter data/bacName
            strcat(bacteria_name[i], bacteriaName);
            
            //add extension
            strcat(bacteria_name[i],".faa");
        }
        
        //close the file - duh
        fclose(input_file);
    }
}

using namespace tidy;

int main(int argc,char * argv[])
{
    time_t t1 = time(NULL);
    
    //Init constants M1, M2 and M
    Init();
    
    //read the bacteria files.
    ReadInputFile(argv[1]);
    
    //Compares bacterias
    CompareAllBacteria();
    
    time_t t2 = time(NULL);
    printf("time elapsed: %ld seconds\n", t2 - t1);
    
    return 0;
}