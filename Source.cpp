/*
#include <iostream>
#include <fstream>
#include <omp.h>
#include <cmath>

#define eps 1e-6

double** ReadMtxMatrix(int& n, int& m) {

    std::ifstream file("tols4000.mtx");
    int num_row, num_col, num_lines;

    while (file.peek() == '%')
        file.ignore(2048, '\n');

    file >> num_row >> num_col >> num_lines;
    n = num_row;
    m = num_col;
    std::cout << std::endl << "The matrix's size is " << num_col << " x " << num_row << std::endl;

    double **Matrix = new double* [num_col]; // ПЕРЕПРОВЕРИТЬ
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

double* CreateSolutionVector(int row, int col, double** Matrix) {
    srand(time(NULL));
    double* SolutionVector = new double[row];
    for (int i = 0; i < row; i++) {
        SolutionVector[i] = rand() % 129 - 64;
    }
    return SolutionVector;
}

double* CreateResultVector(int row, int col, double** A, double* X) {
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
    omp_set_num_threads(8);
    int row, col, N, i, j, k, n;
    double sum1, sum2, diff, norm;
    bool converge = false;
    norm = 1.0;
    double** A = ReadMtxMatrix(row, col);
    convergenceCondition(A, row);
    double* x = CreateSolutionVector(row, col, A);
    double* x_main = new double[row];

    for (i = 0; i < row; i++) {
        x_main[i] = x[i];
    }
    
    double* b = CreateResultVector(row, col, A, x);
    N = row;

    double* x_old = new double[row];

    for (int i = 0; i < row; i++) {
        x_old[i] = 0;
    }

    // Начальное приближение x
    for (i = 0; i < N; i++) {
        x[i] = 0.0;
    }
    double start = omp_get_wtime();
    k = 0;
    // Распараллеливание цикла по k
#pragma omp parallel private(i, j, sum1, sum2, diff)
    
        while (norm > eps) {
            norm = 0.0;
//#pragma omp parallel reduction(+: norm)

            {
                double diff, sum;

//#pragma omp parallel for schedule(dynamic, 50) private(i) shared(x, diff, b, A)
//#pragma omp nowait
                for (i = 0; i < N; i++)
                {
                    diff = b[i];
                    sum = 0.0;

                    for (j = 0; j < i; j++)
                    {
                        sum += A[i][j] * x[j];
                    }

                    for (j = i + 1; j < N; j++)
                    {
                        sum += A[i][j] * x[j];
                    }

                    double new_x = (diff - sum) / A[i][i];
                    norm += (new_x - x[i]) * (new_x - x[i]);
                    x[i] = new_x;
                }
            }
            norm = sqrt(norm);
            std::cout << norm << std::endl;
            k++;
        

            // Синхронизация потоков
//#pragma omp barrier
        
        }
    double finish = omp_get_wtime() - start;

    // Проверка

    double sum = 0.0;

    for (int i = 0; i < row; i++) {
        sum += x_main[i] - x[i];
    }


    std::cout << "I'm Done! Time: " << finish << ". Iterations: " << k << ". Norm = " << norm << std::endl;
    std::cout << "Check. Sum = " << sum << std::endl;
    for (int i = 0; i < 10; i++)
        std::cout << x_main[i] << " " << x[i] << std::endl;
    
}

*/