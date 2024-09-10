/*
#include <iostream>
#include <fstream>
#include <omp.h>

#define epsilon 1e-6

double** ReadMtxMatrix(int& n, int& m) {

    std::ifstream file("tols4000.mtx");
    int num_row, num_col, num_lines;

    while (file.peek() == '%')
        file.ignore(2048, '\n');

    file >> num_row >> num_col >> num_lines;
    n = num_row;
    m = num_col;
    std::cout << std::endl << "The matrix's size is " << num_col << " x " << num_row << std::endl;

    double** Matrix = new double* [num_col]; // ПЕРЕПРОВЕРИТЬ
    for (int l = 0; l < num_row; l++)
        Matrix[l] = new double[num_row];
    //matrix->non_zero = num_lines;
    for (int i = 0; i < num_row; i++) {
        for (int j = 0; j < num_col; j++) {
            Matrix[i][j] = 0;
        }
    }

    for (int l = 0; l < num_lines; l++) {
        double data;
        int row, col;
        file >> row >> col >> data;
        Matrix[row - 1][col - 1] = data;
        Matrix[col - 1][row - 1] = data;
    }

    file.close();
    return Matrix;
}

int* CreateSolutionVector(int row, int col, double** Matrix) {
    srand(time(NULL));
    int* SolutionVector = new int[row];
    for (int i = 0; i < row; i++) {
        SolutionVector[i] = rand() % 129 - 64;
    }
    return SolutionVector;
}

double* CreateResultVector(int row, int col, double** A, int* X) {
    double* ResultVector = new double[row];
    int sum = 0;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            sum += A[i][j] * X[j];
        }
        ResultVector[i] = sum;
        sum = 0;
    }
    return ResultVector;
}

void convergenceCondition(double** A, int N)
{
    srand(time(NULL));
    for (int i = 0; i < N; ++i)
    {
        double sum = 0;
        for (int j = 0; j < N; ++j)
        {
            if (i != j)
                sum += fabs(A[i][j]);
        }

        //(*this)[i][i] = sum + (1.0 / int(rand() % 10 + 1));
        A[i][i] = sum + 0.2;
    }
}

int main() {
    omp_set_dynamic(0);
    omp_set_num_threads(16);
    int row, col;
    double** A = ReadMtxMatrix(row, col);
    convergenceCondition(A, row);
    int* X = CreateSolutionVector(row, col, A);
    double* B = CreateResultVector(row, col, A, X);

    // Разбивка матрицы А на L и U (нижне- и верхнетреугольную матрицы)

    double** L = new double* [row];
    double** U = new double* [row];

    for (int i = 0; i < row; i++) {
        L[i] = new double[col];
        U[i] = new double[col];
    }


    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (i > j) {
                L[i][j] = A[i][j];
            }
            else if (i < j) {
                U[i][j] = A[i][j];
            }
        }
    }

    // Создание компонент для приближений

    double* x_old = new double[row];
    double* x = new double[row];

    double norm, sum1, sum2;

    //Начальное приближение может быть любым, пусть x[i] = B[i]/A[i][i]

    for (int i = 0; i < row; i++) {
        x_old[i] = 0;
    }

    //Начало алгоритма, k - счётчик итераций

    int k = 0;

    double start = omp_get_wtime();

    do {
        k++;
        int i;
        norm = 0.0;
#pragma omp parallel for private(i) shared(x, x_old, B, A)
        for (i = 0; i < row; i++) {
            x[i] = x_old[i];
        }

        for (int i = 0; i < row; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;

            for (int j = 0; j < i; j++) {
                sum1 += L[i][j] * x[j];
            }

            for (int j = i + 1; j < row; j++) {
                sum2 += U[i][j] * x_old[j];
            }

            x[i] = (B[i] - sum1 - sum2) / A[i][i];
        }

        for (int i = 0; i < row; i++) {
            norm += (x[i] - x_old[i]) * (x[i] - x_old[i]);
        }
        norm = sqrt(norm);
        std::cout << norm << std::endl;
        for (int i = 0; i < row; i++) {
            x_old[i] = x[i];
        }
    } while (norm > epsilon);
    double finish = omp_get_wtime() - start;

    // Проверка

    double sum = 0.0;

    for (int i = 0; i < row; i++) {
        sum += X[i] - x[i];
    }

    std::cout << "I'm Done! Time: " << finish << ". Iterations: " << k << ". Norm = " << norm << std::endl;
    std::cout << "Check. Sum = " << sum;
}

*/