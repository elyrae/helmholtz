#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include <string>
#include "linear_matrix.h"

namespace helmholtz {
    using ddFunction = double (*)(const double, const double);

    void seidel_first_boundary(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                               const double err = 1.0E-4, const size_t max_iterations = 100);
    void seidel_third_boundary(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                               const double err = 1.0E-4, const size_t max_iterations = 100);

    void write_matrix(const LinearMatrix& grid, const std::string& file);
}

#endif