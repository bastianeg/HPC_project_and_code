#include <stdlib.h>

double ***
d_malloc_3d(int m, int n, int k) {

    if (m <= 0 || n <= 0 || k <= 0)
        return NULL;

    double ***array3D = malloc(m * sizeof(double **) +
                               m * n * sizeof(double *) +
                               m * n * k * sizeof(double));
    if (array3D == NULL) {
        return NULL;
    }

    for(int i = 0; i < m; i++) {
        array3D[i] = (double **) array3D + m + i * n ;
    }

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            array3D[i][j] = (double *) array3D + m + m * n + i * n * k + j * k;
        }
    }

    return array3D;
}

