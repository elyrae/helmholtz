#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>

using sourceFunction = double (*)(const double, const double, const double);
const int threads = 2;

class LinearMatrix {
    size_t matrix_rows;    
    size_t matrix_columns;
    std::vector<double> v;
public:
    LinearMatrix(const size_t rows, const size_t columns) : matrix_rows(rows), matrix_columns(columns), v(rows*columns) {}

          double& operator()(const size_t i, const size_t j)       { return v[i*matrix_columns + j]; }
    const double& operator()(const size_t i, const size_t j) const { return v[i*matrix_columns + j]; }
    size_t rows()    const { return matrix_rows; }
    size_t size()    const { return rows(); }
    size_t columns() const { return matrix_columns; }

    void fillZero() {
        for (size_t i = 0; i < rows()*columns(); ++i) 
            v[i] = 0.0;   
    }
};

// void print_matrix(const LinearMatrix &matrix) {
//     const int padd = 12;
//     for (size_t i = 0; i < matrix.rows(); ++i) {
//         for (size_t j = 0; j < matrix.columns(); ++j)
//             std::cout << std::setw(padd) << matrix(i, j);
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;
// }

double source(const double x, const double y, const double k) {
    return 2.0*std::sin(M_PI * y) + k*k*(1.0 - x)*x*std::sin(M_PI * y) + M_PI*M_PI*(1 - x)*x*std::sin(M_PI*y);
}

double exact_solution(const size_t i, const size_t j, const double h) {
    const double x = j*h, y = 1.0 - i*h;
    return (1.0 - x)*x*std::sin(M_PI * y);
}

double max_error(const LinearMatrix& res, const double h) {
    double max = 0.0;
    for (size_t i = 1; i < res.rows() - 1;    ++i)
    for (size_t j = 1; j < res.columns() - 1; ++j)
        if (fabs(res(i, j) - exact_solution(i, j, h)) > max)
            max = fabs(res(i, j) - exact_solution(i, j, h));
    return max;
}

void jacobi_iterations_linear(LinearMatrix& res, const sourceFunction F, const double h, const double k, 
                              const double err = 1.0E-4, const size_t max_iterations = 100) {
    LinearMatrix next(res.rows(), res.columns());

    size_t iteration = 0;
    double diff = 0.0;
    do {
        for (size_t i = 1; i < res.rows()    - 1; ++i)
        for (size_t j = 1; j < res.columns() - 1; ++j)
            next(i, j) = ( F(j*h, 1.0 - i*h, k)*h*h + res(i + 1, j) + res(i - 1, j) + res(i, j - 1) + res(i, j + 1) ) / ( 4.0 + h*h*k*k );

        diff = 0.0;
        for (size_t i = 1; i < res.rows()    - 1; ++i)
            for (size_t j = 1; j < res.columns() - 1; ++j) {
                if (fabs(next(i, j) - res(i, j)) > diff)
                    diff = fabs(res(i, j) - next(i, j));

                res(i, j) = next(i, j);
            }

        iteration = iteration + 1;
    } while ( (diff > err) && (iteration < max_iterations) );
}

void seidel_iterations_linear(LinearMatrix& res, const sourceFunction F, const double h, const double k, 
                              const double err = 1.0E-4, const size_t max_iterations = 100) {
    LinearMatrix next(res.rows(), res.columns());

    size_t iteration = 0;
    double diff = 0.0;
    do {
        for (size_t i = 1;           i < res.rows()    - 1; i = i + 1)
        for (size_t j = 1 + (i % 2); j < res.columns() - 1; j = j + 2)
            next(i, j) = ( F(j*h, 1.0 - i*h, k)*h*h + res(i + 1, j) + res(i - 1, j) + res(i, j - 1) + res(i, j + 1) ) / (4.0 + h*h*k*k);

        for (size_t i = 1;                 i < res.rows()    - 1; i = i + 1)
        for (size_t j = 1 + ((i - 1) % 2); j < res.columns() - 1; j = j + 2)
            next(i, j) = ( F(j*h, 1.0 - i*h, k)*h*h + next(i + 1, j) + next(i - 1, j) + next(i, j - 1) + next(i, j + 1) ) / (4.0 + h*h*k*k);

        diff = 0.0;
        for (size_t i = 1; i < res.rows()    - 1; ++i)
            for (size_t j = 1; j < res.columns() - 1; ++j) {
                if (fabs(next(i, j) - res(i, j)) > diff)
                    diff = fabs(res(i, j) - next(i, j));
            
                res(i, j) = next(i, j);
            }

        iteration = iteration + 1;
    } while ( (diff > err) && (iteration < max_iterations) );
}

void seidel_iterations_parallel(LinearMatrix& res, const sourceFunction F, const double h, const double k, 
                                const double err = 1.0E-4, const size_t max_iterations = 100) {
    LinearMatrix next(res.rows(), res.columns());

    size_t iteration = 0;
    double diff = 0.0;
    do {
        #pragma omp parallel for num_threads(threads)
        for (size_t i = 1;           i < res.rows()    - 1; i = i + 1)
        for (size_t j = 1 + (i % 2); j < res.columns() - 1; j = j + 2)
            next(i, j) = ( F(j*h, 1.0 - i*h, k)*h*h + res(i + 1, j) + res(i - 1, j) + res(i, j - 1) + res(i, j + 1) ) / (4.0 + h*h*k*k);

        #pragma omp parallel for num_threads(threads)
        for (size_t i = 1;                 i < res.rows()    - 1; i = i + 1)
        for (size_t j = 1 + ((i - 1) % 2); j < res.columns() - 1; j = j + 2)
            next(i, j) = ( F(j*h, 1.0 - i*h, k)*h*h + next(i + 1, j) + next(i - 1, j) + next(i, j - 1) + next(i, j + 1) ) / (4.0 + h*h*k*k);

        diff = 0.0;
        #pragma omp parallel for num_threads(threads) reduction(max : diff)
        for (size_t i = 1; i < res.rows() - 1; ++i)
            for (size_t j = 1; j < res.columns() - 1; ++j) {
                if (fabs(next(i, j) - res(i, j)) > diff)
                    diff = fabs(res(i, j) - next(i, j));

                res(i, j) = next(i, j);
            }

        iteration = iteration + 1;
    } while ( (diff > err) && (iteration < max_iterations) );
}

void jacobi_iterations_parallel(LinearMatrix& res, const sourceFunction F, const double h, const double k, 
                                const double err = 1.0E-4, const size_t max_iterations = 100) {
    LinearMatrix next(res.rows(), res.columns());

    size_t iteration = 0;
    double diff = 0.0;
    do {
        #pragma omp parallel for num_threads(threads)
        for (size_t i = 1; i < res.rows()    - 1; ++i)
        for (size_t j = 1; j < res.columns() - 1; ++j)
            next(i, j) = ( F(j*h, 1.0 - i*h, k)*h*h + res(i + 1, j) + res(i - 1, j) + res(i, j - 1) + res(i, j + 1) ) / ( 4.0 + h*h*k*k );

        diff = 0.0;
        #pragma omp parallel for num_threads(threads) reduction(max : diff)
        for (size_t i = 1; i < res.rows()    - 1; ++i)
            for (size_t j = 1; j < res.columns() - 1; ++j) {
                if (fabs(next(i, j) - res(i, j)) > diff)
                    diff = fabs(res(i, j) - next(i, j));

                res(i, j) = next(i, j);
            }

        iteration = iteration + 1;
    } while ( (diff > err) && (iteration < max_iterations) );
}

int main() {
    const size_t N = 100;
    const double h = 1.0 / double(N - 1);
    const double k = 2.0;
    LinearMatrix result(N, N);

    double start_sequential = omp_get_wtime();
    seidel_iterations_linear(result, source, h, k, 1.0E-7, 20000);
    //jacobi_iterations_linear(result, source, h, k, 1.0E-7, 20000);
    double end_sequential   = omp_get_wtime();
    std::cout << "max error sequential: " << max_error(result, h) << "\n";

    result.fillZero();

    double start_parallel = omp_get_wtime();
    seidel_iterations_parallel(result, source, h, k, 1.0E-7, 10000);
    //jacobi_iterations_parallel(result, source, h, k, 1.0E-7, 20000);
    double end_parallel   = omp_get_wtime();
    std::cout << "max error parallel:   " << max_error(result, h) << "\n";

    const double acceleration = (end_sequential - start_sequential) / (end_parallel - start_parallel);
    std::cout << std::endl;
    std::cout << "sequential algorithm:      " << (end_sequential - start_sequential) << std::endl;
    std::cout << "parallel algorithm:        " << (end_parallel - start_parallel) << std::endl;
    std::cout << "acceleration (efficiency): " << acceleration << " (" << (acceleration / threads) << ")" << std::endl;
    return 0;
}
