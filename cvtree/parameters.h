/*
 A helper header that contains structs used to pass parameters to threads.
 
 Luiz Ribeiro (n9383298) QUT Parallel Computing // 2015 - 2nd semester
*/

#ifndef cvtree_parameters_h
#define cvtree_parameters_h

struct readingParams {
    char *file_name;
    int  thread_id;
};

#endif
