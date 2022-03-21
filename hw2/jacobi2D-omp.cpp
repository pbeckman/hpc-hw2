// Poisson 2D solver using Jacobi
// $ g++ -std=c++11 -O3 -o jacobi2D-omp jacobi2D-omp.cpp && ./jacobi2D-omp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

void solve(int N, int maxiter, double* f, double* u, double* u_l) {
  double h2 = pow(1.0/(N+1), 2);

  for (long iter = 0; iter < maxiter; iter++) {
    #pragma omp parallel for collapse(2)
    for (long j = 1; j < N+1; j++) {
      for (long i = 1; i < N+1; i++) {
        u[i+j*N] = (h2*f[i+j*N] + u_l[i-1+j*N] + u_l[i+1+j*N] + u_l[i+(j-1)*N] + u_l[i+(j+1)*N]) / 4;
      }
    }

    #pragma omp parallel for collapse(2)
    for (long j = 1; j < N+1; j++) {
      for (long i = 1; i < N+1; i++) {
        u_l[i+j*N] = u[i+j*N];
      }
    }
  }
}

int main(int argc, char** argv) {
  Timer t;
  int maxiter = 100;

  printf("         N       Time\n");
  for (int N = 8; N <= 8192; N *= 2) {
    // allocate necessary arrays
    double* f   = (double*) malloc((N+2)*(N+2)*sizeof(double));
    double* u   = (double*) malloc((N+2)*(N+2)*sizeof(double));
    double* u_l = (double*) malloc((N+2)*(N+2)*sizeof(double));

    // initialize arrays
    for (long j = 1; j < N+1; j++) {
      for (long i = 1; i < N+1; i++) {
        f[i+j*N]   = 1;
        u[i+j*N]   = 0; 
        u_l[i+j*N] = 0;
      }
    }

    // perform maxiter iterations of solver
    t.tic();
    solve(N, maxiter, f, u, u_l);
    printf("%10i %10f\n", N, t.toc());
    
    free(f);
    free(u);
    free(u_l);
  }
}