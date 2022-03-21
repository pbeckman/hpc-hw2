// Poisson 2D solver using Gauss-Seidel
// $ g++ -std=c++11 -O3 -o gs2D-omp gs2D-omp.cpp && ./gs2D-omp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

void solve(int N, int maxiter, double* f, double* u) {
  double h2 = pow(1.0/(N+1), 2);

  for (long iter = 0; iter < maxiter; iter++) {
    #pragma omp parallel for
    // red sweep
    for (long j = 1; j < N+1; j++) {
      for (long i = j%2; i < N+1; i += 2) {
        u[i+j*N] = (h2*f[i+j*N] + u[i-1+j*N] + u[i+1+j*N] + u[i+(j-1)*N] + u[i+(j+1)*N]) / 4;
      }
    }

    #pragma omp parallel for
    // black sweep
    for (long j = 1; j < N+1; j++) {
      for (long i = 1+(j%2); i < N+1; i += 2) {
        u[i+j*N] = (h2*f[i+j*N] + u[i-1+j*N] + u[i+1+j*N] + u[i+(j-1)*N] + u[i+(j+1)*N]) / 4;
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

    // initialize arrays
    for (long j = 1; j < N+1; j++) {
      for (long i = 1; i < N+1; i++) {
        f[i+j*N]   = 1;
        u[i+j*N]   = 0; 
      }
    }

    // perform maxiter iterations of solver
    t.tic();
    solve(N, maxiter, f, u);
    printf("%10i %10f\n", N, t.toc());
    
    free(f);
    free(u);
  }
}