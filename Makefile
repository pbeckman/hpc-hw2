all: val_test mmult omp_bug # poisson

val_test:
	g++ -O0 -g -march=native -o val_test01_solved val_test01_solved.cpp
	g++ -O0 -g -march=native -o val_test02_solved val_test02_solved.cpp

mmult: 
	g++ -std=c++11 -fopenmp -O3 -march=native -o MMult1 MMult1.cpp
	
omp_bug:
	g++ -std=c++11 -fopenmp -O3 -march=native -o omp_solved2 omp_solved2.c
	g++ -std=c++11 -fopenmp -O3 -march=native -o omp_solved3 omp_solved3.c
	g++ -std=c++11 -fopenmp -O3 -march=native -o omp_solved4 omp_solved4.c
	g++ -std=c++11 -fopenmp -O3 -march=native -o omp_solved5 omp_solved5.c
	g++ -std=c++11 -fopenmp -O3 -march=native -o omp_solved6 omp_solved6.c

# poisson:
# 	g++ -std=c++11 -fopenmp -O3 -march=native -o jacobi2D-omp jacobi2D-omp.cpp
# 	g++ -std=c++11 -fopenmp -O3 -march=native -o gs2D-omp gs2D-omp.cpp