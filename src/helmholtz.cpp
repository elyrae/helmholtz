#include "helmholtz.h"

#include <iostream>
#include <fstream>
#include <cmath>

void helmholtz::writeMatrix(const LinearMatrix& grid, const double h, const std::string& filename)
{
    std::ofstream out(filename);

    double x = 0.0, y = 0.0;
    for (size_t i = 0; i < grid.rows()   ; ++i)
    for (size_t j = 0; j < grid.columns(); ++j) {
        x = j*h;
        y = 1.0 - i*h;

        out << x << " " << y << " " << grid(i, j) << "\n";
    }
}

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

        // std::cout << iteration << ": " << diff << std::endl;
        if (iteration < 10) {
            helmholtz::writeMatrix(next, h, std::to_string(iteration + 1) + ".txt");
        }

        iteration++;
    } while ( (diff > err) && (iteration < max_iterations) );
}

void helmholtz::seidel(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                       const double err, const size_t max_iterations) 
{
    LinearMatrix next(m.rows(), m.columns());

    const double h = 1.0 / double(m.size() - 1);

    size_t iteration = 0;
    double diff = 0.0;
    double x = 0.0, y = 0.0, a_x_forw = 0.0, a_x_back = 0.0, a_y_forw = 0.0, a_y_back = 0.0;    
    do {
        for (size_t i = 1;           i < m.rows()    - 1; i = i + 1)
        for (size_t j = 1 + (i % 2); j < m.columns() - 1; j = j + 2) {
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

        for (size_t i = 1;                 i < m.rows()    - 1; i = i + 1)
        for (size_t j = 1 + ((i - 1) % 2); j < m.columns() - 1; j = j + 2) {
            x = j*h;
            y = 1.0 - i*h;

            a_x_forw = lamb(x + h/2.0, y);
            a_x_back = lamb(x - h/2.0, y);
            a_y_forw = lamb(x, y + h/2.0);
            a_y_back = lamb(x, y - h/2.0);

            next(i, j) = (Q(x, y)*h*h + next(i + 1, j)*a_x_forw + next(i - 1, j)*a_x_back 
                                      + next(i, j + 1)*a_y_forw + next(i, j - 1)*a_y_back) 
                       / (a_x_forw + a_x_back + a_y_forw + a_y_back + k(x, y)*h*h);
        }

        diff = 0.0;
        for (size_t i = 1; i < m.rows()    - 1; ++i)
        for (size_t j = 1; j < m.columns() - 1; ++j) {
            if (fabs(next(i, j) - m(i, j)) > diff)
                diff = fabs(m(i, j) - next(i, j));
            m(i, j) = next(i, j);
        }

        std::cout << diff << std::endl;

        iteration++;
    } while ( (diff > err) && (iteration < max_iterations) );
}
