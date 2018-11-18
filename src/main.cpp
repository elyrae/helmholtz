#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>

#include "helmholtz.h"
#include "linear_matrix.h"

double lamb(const double x, const double)
{
    // return 1.0 + x;
    return 1.0;
}

double k(const double, const double)
{
    return 2.0;
}

double Q(const double x, const double y) {
    return -2*pow(-1 + y,2)*pow(y,2) + 12*x*pow(-1 + y,2)*pow(y,2) - 4*pow(x,3)*(-1 + 6*y - 5*pow(y,2) - 2*pow(y,3) + pow(y,4)) + 
   2*pow(x,4)*(-1 + 6*y - 5*pow(y,2) - 2*pow(y,3) + pow(y,4)) - 2*pow(x,2)*(1 - 6*y + 11*pow(y,2) - 10*pow(y,3) + 5*pow(y,4));

    // return -2*pow(-1 + y,2)*pow(y,2) + 12*x*pow(-1 + y,2)*pow(y,2) - 4*pow(x,3)*(-1 + 6*y - 5*pow(y,2) - 2*pow(y,3) + pow(y,4)) + 
    // 2*pow(x,4)*(-1 + 6*y - 5*pow(y,2) - 2*pow(y,3) + pow(y,4)) - 2*pow(x,2)*(1 - 6*y + 11*pow(y,2) - 10*pow(y,3) + 5*pow(y,4));   

    // return 10.0;
    //return -1.0 - 10.0*x + 4.0*x*x - 4.0*y + 4.0*y*y;
    // return -2.0 - 4.0*x + 4.0*x*x - 4.0*y + 4.0*y*y;
    // return 2.0*std::sin(M_PI * y) + k(x, y)*(1.0 - x)*x*std::sin(M_PI * y) + M_PI*M_PI*(1 - x)*x*std::sin(M_PI*y);
}

// double exact_solution(const size_t i, const size_t j, const double h) {
//     const double x = j*h, y = 1.0 - i*h;
//     return (1.0 - x)*x*std::sin(M_PI * y);
// }

// double max_error(const LinearMatrix& m) {
//     const double h = 1.0 / double(m.size() - 1);
//     double max = 0.0;
//     for (size_t i = 1; i < m.rows()    - 1; ++i)
//     for (size_t j = 1; j < m.columns() - 1; ++j)
//         if (fabs(m(i, j) - exact_solution(i, j, h)) > max)
//             max = fabs(m(i, j) - exact_solution(i, j, h));
//     return max;
// }

int main() {
    const size_t N = 160;
    // const double h = 1.0 / double(N - 1);
    LinearMatrix m(N, N);

    // double x = 0.0, y = 0.0;
    for (size_t i = 0; i < m.rows(); ++i) {
        // x = i*h;
        // y = i*h;

        m(0,            i) = 0.0; // 0.5 - y + y*y;
        m(m.size() - 1, i) = 0.0; // 0.5 - y + y*y;
        m(i,            0) = 0.0; // 0.5 - x + x*x;
        m(i, m.size() - 1) = 0.0; // 0.5 - x + x*x; 
    }

    // helmholtz::writeMatrix(m, "0.txt");
    // helmholtz::jacobi(m, lamb, k, Q, 1.0E-4, 1);

    // helmholtz::jacobi(m, lamb, k, Q, 1.0E-4, 1000);
    helmholtz::jacobiThirdBoundary(m, lamb, k, Q, 1.0E-4, 50000);

    return 0;
}
