#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

double pi_int(int N){
    int i;
    double Nr = 1.0/N;
    double pi = 0.0;
    #pragma omp parallell for default(none) \
            shared(Nr,N) private(i,pi) \
            reduction(+: pi)
    for(i=1;i<=N;i++){
        pi += 4/(1+Nr*Nr*(i-0.5)*(i-0.5));

    }

    pi *= Nr;
    
    return pi;
}

int main(int argc, char *argv[]){
    int N = 1000000000;
    
    double pi = pi_int(N);
    printf("Pi is %.5f with %d iterations\n",pi,N);
    
    return 0;
}
