#include <stdio.h>

double pi_int(int N){
    int i;
    double Nr = 1.0/N;
    double pi = 0.0;
    double tmp;
    for(i=0;i<N;i++){
        tmp = 4/(1+Nr*Nr*(i-0.5)*(i-0.5));
        pi += tmp;
    }
    pi *= Nr;
    return pi;
}

int main(int argc, char *argv[]){
    int N = 10000;
    double pi = pi_int(N);
    printf("Pi is %.5f with %d iterations\n",pi,N);
    
    return 0;
}
