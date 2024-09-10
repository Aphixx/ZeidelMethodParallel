#include <omp.h> 
#include <iostream> 

int test() {
	float a = omp_get_wtime();
	omp_get_thread_num();
	omp_set_dynamic(0);
	omp_set_num_threads(5);
	int sum = 0;
	int* A = new int[1000000000];
	for (int i = 0; i < 1000000; i++)
	{
		*(A + i) = 1;
	}
	double start = omp_get_wtime();
#pragma omp parallel for 
	for (int i = 0; i < 1000000; i++)
	{
#pragma omp critical 
		sum += *(A + i);
	}
	double finish = omp_get_wtime() - start;
	std::cout << sum << std::endl;
	std::cout << a << std::endl;
	std::cout << finish;
	return 0;
}