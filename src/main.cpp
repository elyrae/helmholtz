#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>

#include "helmholtz.h"
#include "linear_matrix.h"

double lambda(const double x, const double)
{
    // return 1.0 + x;
    return 1.0;
}

double kappa(const double, const double)
{
    return 2.0;
}

double source(const double x, const double y) {
    return (-2*(2399 - 5000*x + 5000*pow(x,2) - 5000*y + 5000*pow(y,2)))/pow(M_E,50*(pow(-0.5 + x,2) + pow(-0.5 + y,2)));

    // return -2*pow(-1 + y,2)*pow(y,2) + 12*x*pow(-1 + y,2)*pow(y,2) - 4*pow(x,3)*(-1 + 6*y - 5*pow(y,2) - 2*pow(y,3) + pow(y,4)) + 
    // 2*pow(x,4)*(-1 + 6*y - 5*pow(y,2) - 2*pow(y,3) + pow(y,4)) - 2*pow(x,2)*(1 - 6*y + 11*pow(y,2) - 10*pow(y,3) + 5*pow(y,4));

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
    const size_t NX = 50;
    const size_t NY = 100;
    const double hx = 1.0 / double(NX - 1);
    const double hy = 1.0 / double(NY - 1);
    LinearMatrix m(NX, NY), k(NX, NY), lamb(NX, NY), Q(NX, NY);

    double x = 0.0, y = 0.0;
    for (size_t i = 0; i < m.rows();    ++i)
    for (size_t j = 0; j < m.columns(); ++j) {
        x = i*hx;
        y = j*hy;

        k(i, j) = kappa(x, y);
        lamb(i, j) = lambda(x, y);
        Q(i, j) = source(x, y);
    }

    for (size_t i = 0; i < m.rows(); ++i) {
        // x = i*h;
        // y = i*h;

        m(0,            i) = 0.0; // 0.5 - y + y*y;
        m(m.size() - 1, i) = 0.0; // 0.5 - y + y*y;
        m(i,            0) = 0.0; // 0.5 - x + x*x;
        m(i, m.size() - 1) = 0.0; // 0.5 - x + x*x; 
    }

    double start = omp_get_wtime();
    helmholtz::seidel_third_boundary(m, lamb, k, Q, 1.0E-4, 1000);
    double end = omp_get_wtime();

    printf("Time: %f s\n", end - start);

    return 0;
}
