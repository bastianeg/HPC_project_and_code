#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char *argv[]) {
    int t_id = 0;
    #pragma omp parallel private(t_id) num_threads(4)
    {
        t_id = omp_get_thread_num();
        printf("Hello world from %d!\n",t_id);
    }
    return(0);
}

