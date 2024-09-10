#include <iostream>
#include <omp.h>


// omp_get_wtime()
// omp_get_num_threads()
// omp_get_thread_num()
// omp_set_dynamic(0)
// omp_set_num_threads(5)


double multijk(float** result, float** left, float** right, int m, int n, int p) {
	double start = omp_get_wtime();
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < p; ++j) {
			for (int k = 0; k < n; ++k) {
				result[i][j] += left[i][k] * right[k][j];
			}
		}
	}
	return omp_get_wtime() - start;
}

double multikj(float** result, float** left, float** right, int m, int n, int p) {
	double start = omp_get_wtime();
	for (int i = 0; i < m; ++i) {
		for (int k = 0; k < n; ++k) {
			for (int j = 0; j < p; ++j) {
				result[i][j] += left[i][k] * right[k][j];
			}
		}
	}
	return omp_get_wtime() - start;
}

double multjik(float** result, float** left, float** right, int m, int n, int p) {
	double start = omp_get_wtime();
	for (int j = 0; j < p; ++j) {
		for (int i = 0; i < m; ++i) {
			for (int k = 0; k < n; ++k) {
				result[i][j] += left[i][k] * right[k][j];
			}
		}
	}
	return omp_get_wtime() - start;
}

double multjki(float** result, float** left, float** right, int m, int n, int p) {
	double start = omp_get_wtime();
	for (int j = 0; j < p; ++j) {
		for (int k = 0; k < n; ++k) {
			for (int i = 0; i < m; ++i) {
				result[i][j] += left[i][k] * right[k][j];
			}
		}
	}
	return omp_get_wtime() - start;
}

double multkij(float** result, float** left, float** right, int m, int n, int p) {
	double start = omp_get_wtime();
	for (int k = 0; k < n; ++k) {
		for (int i = 0; i < m; ++i) {
#pragma omp for
			for (int j = 0; j < p; ++j)
				result[i][j] += left[i][k] * right[k][j];
		}
	}
	return omp_get_wtime() - start;
}

double multkji(float** result, float** left, float** right, int m, int n, int p) {
	double start = omp_get_wtime();
	for (int k = 0; k < n; ++k) {
		for (int j = 0; j < p; ++j) {
			for (int i = 0; i < m; ++i) {
				result[i][j] += left[i][k] * right[k][j];
			}
		}
	}
	return omp_get_wtime() - start;
}

void print(float** matrix, int m, int n) {
	std::cout << "\n" << "Matrix:\n";
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

float** gen_matrix(int m, int n, float fill) {
	float** matrix = new float* [m];
	for (int i = 0; i < m; ++i) {
		matrix[i] = new float[n];
		for (int j = 0; j < n; ++j) {
			matrix[i][j] = fill;
		}
	}
	return matrix;
}

int lection() {
	int m = 0;
	int n = 0;
	int p = 0;

	float** left = gen_matrix(1000, 2000, 1);

	float** right = gen_matrix(2000, 1000, 2);

	float** result = gen_matrix(1000, 1000, 0);

	std::cout << multkij(result, left, right, 1000, 2000, 1000);

	print(result, 1, 10);

	return 0;
}