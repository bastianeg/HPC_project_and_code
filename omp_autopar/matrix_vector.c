#include <stdio.h>
#include <stdlib.h>
void mvx(int m, int n, double *a, double *b, double *c){
    int i,j;
    double sum;

    #pragma omp parallel for default(none) \
    shared(m,n,a,b,c) private(i,j,sum)
    for(i=0;i<m;i++){
        sum = 0.0;
        for(j=0;j<n;j++){
            sum += b[i*n+j] * c[j];
        }
        a[i] = sum;
    }
}

int main(int argc, char* argv[]){
    int i,j;
    int n = 1000;
    int m = n;
    double *a = malloc(m*sizeof(double));
    double *c = malloc(n*sizeof(double));
    double *b = malloc(n*m*sizeof(double));
    
    for(i=0;i<m;i++){    
        for(j=0;j<n;j++){
            b[i*n+j] = i*j;
        }
    }
    for(j=0;j<n;j++){
        c[j] = j;
    }
    mvx(m,n,a,b,c);
    
    for(i=0;i<m;i++){
        printf("%.1f\n",a[i]);
    }

    free(a);
    free(b);
    free(c);
    return 0;
}