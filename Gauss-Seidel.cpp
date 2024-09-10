#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;

double** ReadMtxMatrix(int& n, int& m) {

    std::ifstream file("tols1090.mtx");
    int num_row, num_col, num_lines;

    while (file.peek() == '%')
        file.ignore(2048, '\n');

    file >> num_row >> num_col >> num_lines;
    n = num_row;
    m = num_col;
    std::cout << std::endl << "The matrix's size is " << num_col << " x " << num_row << std::endl;

    double** Matrix = new double* [num_col]; // оепеопнбепхрэ
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

void writeBVector(double* b, int N)
{
    std::ofstream file("B.txt");

    file << N << "\n";

    for (int i = 0; i < N; ++i)
        file << b[i] << " ";
}

double* readTXTVector()
{
    std::ifstream file("B.txt");
    int size;

    file >> size;
    double* b = new double[size];

    for (int i = 0; i < size; ++i)
        file >> b[i];

    return b;
}

int main() {
    const double tolerance = 1e-9; // convergence tolerance

    double fulltime = omp_get_wtime();

    int row, col, N, i, j, k, n;
    double sum1, sum2, diff, norm;
    bool converge = false;
    norm = 1.0;
    double** A = ReadMtxMatrix(row, col);
    convergenceCondition(A, row);
    double* x = CreateSolutionVector(row, col, A);
    //double* b = CreateResultVector(row, col, A, x);
    //writeBVector(b, row);
    double* b = readTXTVector();

    //cout << b[0] << " " << b[row] << endl;

    double* x_main = new double[row];

    for (i = 0; i < row; i++) {
        x_main[i] = x[i];
    }

    double sr = 0;

    for (i = 0; i < row; i++) {
        for (j = 0; j < row; j++) {
            sr += A[i][j]*A[i][j];
        }
    }

    //cout << sqrt(sr) / 8784 << endl;

    N = row;
    double* x_new = new double[row];
    for (i = 0; i < N; i++) {
        x[i] = 0.0;
    }
    // initialize A, b, and x

    // set number of threads to use
    int num_threads = omp_get_max_threads();
    //cout << num_threads;
    omp_set_dynamic(0);
    omp_set_num_threads(8);

    
    double start = omp_get_wtime();
    bool converged = false;
    int iter = 0;
    while (!converged) {
#pragma omp parallel
        {
#pragma omp for
            for (int i = 0; i < N; i++) {
                //cout << omp_get_num_threads() << endl;
                double sum = 0.0;
                for (int j = 0; j < N; j++) {
                    if (j != i) {
                        sum += A[i][j] * x[j];
                    }
                }
                x_new[i] = (b[i] - sum) / A[i][i];
            }

#pragma omp barrier

#pragma omp single
            {
                double error = 0.0;
                for (int i = 0; i < N; i++) {
                    double diff = abs(x_new[i] - x[i]);
                    //cout << diff << endl;
                    if (diff > error) {
                        error = diff;
                    }
                    x[i] = x_new[i];
                }
                iter++;
                if (error < tolerance) {
                    converged = true;
                }
            }
        }
    }
    double finish = omp_get_wtime();
    // output solution and number of iterations
   

    std::cout << "I'm Done! Time: " << finish - start << ". Iterations: " << iter << std::endl;
    
    
    /*
    for (int i = 0; i < 7; i++) {
        cout << x_main[i] << " " << x[i] << endl;
    }
    */
    return 0;
}
