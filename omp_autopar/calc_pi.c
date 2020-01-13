#include <stdio.h>

int main(int argc, char *argv[]){
    int i;
    int N = 10000;
    double Nr = 1.0/N;
    double pi = 0.0;
    double tmp;
    for(i=0;i<N;i++){
        tmp = 4/(1+Nr*Nr*(i-0.5)*(i-0.5));
        pi += tmp;
    }
    pi *= Nr;
    printf("Pi is %.5f with %d iterations\n",pi,N);
    
    return 0;
}