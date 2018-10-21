#include "helmholtz.h"

#include <iostream>
#include <cmath>

// double source(const double x, const double y, const double k) {
//     return 2.0*std::sin(M_PI * y) + k*k*(1.0 - x)*x*std::sin(M_PI * y) + M_PI*M_PI*(1 - x)*x*std::sin(M_PI*y);
// }

// double exact_solution(const size_t i, const size_t j, const double h) {
//     const double x = j*h, y = 1.0 - i*h;
//     return (1.0 - x)*x*std::sin(M_PI * y);
// }

// double max_error(const LinearMatrix& res, const double h) {
//     double max = 0.0;
//     for (size_t i = 1; i < res.rows() - 1;    ++i)
//     for (size_t j = 1; j < res.columns() - 1; ++j)
//         if (fabs(res(i, j) - exact_solution(i, j, h)) > max)
//             max = fabs(res(i, j) - exact_solution(i, j, h));
//     return max; λ
// }

void helmholtz::jacobi(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                       const double err, const size_t max_iterations) 
{
    LinearMatrix next(m.rows(), m.columns());

    const double h = 1.0 / double(m.size() - 1);

    size_t iteration = 0;
    double diff = 0.0;
    double x = 0.0, y = 0.0, a_x_forw = 0.0, a_x_back = 0.0, a_y_forw = 0.0, a_y_back = 0.0;
    do {
        for (size_t i = 1; i < m.rows()    - 1; ++i)
        for (size_t j = 1; j < m.columns() - 1; ++j) {
            x = j*h;
            y = 1.0 - i*h;

            a_x_forw = lamb(x + h/2.0, y);
            a_x_back = lamb(x - h/2.0, y);
            a_y_forw = lamb(x, y + h/2.0);
            a_y_back = lamb(x, y - h/2.0);

            next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back 
                                      + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
                       / (a_x_forw + a_x_back + a_y_forw + a_y_back + k(x, y)*h*h);
        }

        diff = 0.0;
        for (size_t i = 1; i < m.rows()    - 1; ++i)
        for (size_t j = 1; j < m.columns() - 1; ++j) {
           if (fabs(next(i, j) - m(i, j)) > diff)
               diff = fabs(next(i, j) - m(i, j));
           m(i, j) = next(i, j);
        }

        std::cout << diff << std::endl;

        iteration++;
    } while ( (diff > err) && (iteration < max_iterations) );
}
