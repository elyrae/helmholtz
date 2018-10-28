#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include "linear_matrix.h"

#include <string>

namespace helmholtz {
    using ddFunction = double (*)(const double, const double);

    void jacobi(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                const double err = 1.0E-4, const size_t max_iterations = 100);
    void seidel(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                const double err = 1.0E-4, const size_t max_iterations = 100);

    void writeMatrix(const LinearMatrix& grid, const double h, const std::string& file);
}

#endif