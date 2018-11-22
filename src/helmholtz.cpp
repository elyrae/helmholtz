#include "helmholtz.h"

#include <iostream>
#include <fstream>
#include <cmath>

const int THREADS = 4;

void helmholtz::write_matrix(const LinearMatrix& grid, const std::string& filename)
{
    std::ofstream out(filename);

    for (size_t i = 0; i < grid.rows(); ++i) {
        for (size_t j = 0; j < grid.columns(); ++j)
            out << grid(i, j) << " ";
        out << "\n"; 
    }
}

void helmholtz::seidel_first_boundary(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                                      const double err, const size_t max_iterations) 
{
    LinearMatrix next(m.rows(), m.columns());
    for (size_t i = 0; i < m.rows(); ++i) {
        next(0,            i) = m(0,            i);
        next(m.size() - 1, i) = m(m.size() - 1, i);
        next(i,            0) = m(i,            0);
        next(i, m.size() - 1) = m(i, m.size() - 1); 
    }

    const double h = 1.0 / double(m.size() - 1);

    size_t iteration = 0;
    double x = 0.0, y = 0.0, a_x_forw = 0.0, a_x_back = 0.0, a_y_forw = 0.0, a_y_back = 0.0;
    do {
        for (size_t i = 1; i < m.rows()    - 1; ++i)
        for (size_t j = 1; j < m.columns() - 1; ++j) {
            x = i*h;
            y = j*h;

            a_x_forw = lamb(x + h/2.0, y);
            a_x_back = lamb(x - h/2.0, y);
            a_y_forw = lamb(x, y + h/2.0);
            a_y_back = lamb(x, y - h/2.0);
            next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back 
                                      + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
                         / (a_x_forw + a_x_back + a_y_forw + a_y_back + k(x, y)*h*h);
        }

        // diff = 0.0;
        // for (size_t i = 1; i < m.rows()    - 1; ++i)
        // for (size_t j = 1; j < m.columns() - 1; ++j) {
            // if (fabs(next(i, j) - m(i, j)) > diff)
            //    diff = fabs(next(i, j) - m(i, j));        }

        // std::cout << iteration << ": " << diff << std::endl;
        // if ( ((iteration + 1) % 10) == 0 ) {
        //     helmholtz::writeMatrix(next, std::to_string(iteration + 1) + ".txt");
        // }

        m.swap(next);
        iteration++;
        if ( (iteration % 100) == 0 ) {
            helmholtz::write_matrix(next, std::to_string(iteration) + ".txt");
        }        
    } while ( /*(diff > err) &&*/ (iteration < max_iterations) );
}

double exact(const double x, const double y)
{
    return exp(-50.0 * ( (x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5) ));
}

void seidel_step(const size_t i, const size_t j, LinearMatrix& next, LinearMatrix& m, 
                 const helmholtz::ddFunction lamb, const helmholtz::ddFunction k, const helmholtz::ddFunction Q)
{
    const double hx = 1.0 / double(m.rows()    - 1);
    const double hy = 1.0 / double(m.columns() - 1);
    const double hy_hx = hy / hx;
    const double hx_hy = hx / hy;

    const double x = i*hx;
    const double y = j*hy;

    const double a_x_forw = (i != (m.rows() - 1)   )*lamb(x + hx/2.0, y         );
    const double a_x_back = (i != 0                )*lamb(x - hx/2.0, y         );
    const double a_y_forw = (j != (m.columns() - 1))*lamb(x,          y + hy/2.0);
    const double a_y_back = (j != 0                )*lamb(x,          y - hy/2.0);

    // =========================================================== 
    // внутренняя часть области
    if ((0 < i) && (i < (m.rows() - 1)) && (0 < j) && (j < (m.columns() - 1)))
        next(i, j) = (Q(x, y)*hy*hy + m(i + 1, j)*a_x_forw*hy_hx*hy_hx 
                                    + m(i - 1, j)*a_x_back*hy_hx*hy_hx 
                                    + m(i, j + 1)*a_y_forw 
                                    + m(i, j - 1)*a_y_back) 
                      / ((a_x_forw + a_x_back)*hy_hx*hy_hx + a_y_forw + a_y_back + k(x, y)*hy*hy);

    // ===========================================================
    // левая граница без угловых точек
    if ((i == 0) && (0 < j) && (j < (m.columns() - 1)))
        next(i, j) = (Q(x, y)*hy*hy*0.5 + m(1, j    )*a_x_forw*hy_hx*hy_hx 
                                        + m(0, j + 1)*a_y_forw*0.5 
                                        + m(0, j - 1)*a_y_back*0.5) 
                     / (a_x_forw*hy_hx*hy_hx + 1.5*k(x, y)*lamb(x, y)*hy*hy_hx + 0.5*a_y_forw + 0.5*a_y_back + k(x, y)*hy*hy*0.5);

    // правая граница без угловых точек
    if ((i == (m.rows() - 1)) && (0 < j) && (j < (m.columns() - 1)))
        next(i, j) = (Q(x, y)*hy*hy*0.5 + m(i - 1, j)*a_x_back*hy_hx*hy_hx 
                                        + m(i, j + 1)*a_y_forw*0.5 
                                        + m(i, j - 1)*a_y_back*0.5) 
                     / (a_x_back*hy_hx*hy_hx + 1.5*k(x, y)*lamb(x, y)*hy*hy_hx + 0.5*a_y_forw + 0.5*a_y_back + k(x, y)*hy*hy*0.5);

    // верхняя граница без угловых точек
    if ((0 < i) && (i < (m.rows() - 1)) && (j == (m.columns() - 1)))
        next(i, j) = (Q(x, y)*hx*hx*0.5 + m(i + 1, j)*a_x_forw*0.5 
                                        + m(i - 1, j)*a_x_back*0.5 
                                        + m(i, j - 1)*a_y_back*hx_hy*hx_hy) 
                     / (0.5*a_x_forw + 0.5*a_x_back + 1.5*k(x, y)*lamb(x, y)*hx*hx_hy + a_y_back*hx_hy*hx_hy + k(x, y)*hx*hx*0.5);

    // нижняя граница без угловых точек
    if ((0 < i) && (i < (m.rows() - 1)) && (j == 0))
        next(i, j) = (Q(x, y)*hx*hx*0.5 + m(i + 1, 0)*a_x_forw*0.5 
                                        + m(i - 1, 0)*a_x_back*0.5 
                                        + m(i,     1)*a_y_forw*hx_hy*hx_hy) 
                     / (0.5*a_x_forw + 0.5*a_x_back + a_y_forw*hx_hy*hx_hy + 1.5*k(x, y)*lamb(x, y)*hx*hx_hy + k(x, y)*hx*hx*0.5);

    // ===========================================================                  
    // левый верхний угол +
    if ((i == 0) && (j == (m.columns() - 1)))
        next(i, j) = (Q(x, y)*hy*hy*0.5 + m(1, j    )*a_x_forw*hy_hx*hy_hx 
                                        + m(0, j - 1)*a_y_back)
                     / (a_x_forw*hy_hx*hy_hx + 1.5*k(x, y)*lamb(x, y)*hy + 1.5*k(x, y)*lamb(x, y)*hy*hy_hx + a_y_back + k(x, y)*hy*hy*0.5);

    // левый нижний угол +
    if ((i == 0) && (j == 0))
        next(i, j) = (Q(x, y)*hy*hy*0.5 + m(1, 0)*a_x_forw*hy_hx*hy_hx  
                                        + m(0, 1)*a_y_forw) 
                     / (a_x_forw*hy_hx*hy_hx + 1.5*k(x, y)*lamb(x, y)*hy + a_y_forw + 1.5*k(x, y)*lamb(x, y)*hy*hy_hx + k(x, y)*hy*hy*0.5);

    // правый верхний угол +++
    if ((i == (m.rows() - 1)) && (j == (m.columns() - 1)))
        next(i, j) = (Q(x, y)*hy*hy*0.5 + m(i - 1, j)*a_x_back*hy_hx*hy_hx 
                                        + m(i, j - 1)*a_y_back) 
                     / (1.5*k(x, y)*lamb(x, y)*hy + a_x_back*hy_hx*hy_hx + 1.5*k(x, y)*lamb(x, y)*hy*hy_hx + a_y_back + k(x, y)*hy*hy*0.5); 

    // правый нижний угол
    if ((i == (m.rows() - 1)) && (j == 0)) 
        next(i, j) = (Q(x, y)*hy*hy*0.5 + m(i - 1, 0)*a_x_back*hy_hx*hy_hx 
                                        + m(i,     1)*a_y_forw) 
                     / (1.5*k(x, y)*lamb(x, y)*hy + a_x_back*hy_hx*hy_hx + a_y_forw + 1.5*k(x, y)*lamb(x, y)*hy*hy_hx + k(x, y)*hy*hy*0.5);
}

void helmholtz::seidel_third_boundary(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                                      const double err, const size_t max_iterations) 
{
    LinearMatrix next(m.rows(), m.columns());

    size_t iteration = 0;
    #pragma omp parallel num_threads(THREADS)
    {
        do {
            #pragma omp for
            for (size_t i = 1;           i < m.rows()    - 1; i = i + 1)
            for (size_t j = 1 + (i % 2); j < m.columns() - 1; j = j + 2)
                seidel_step(i, j, next, m, lamb, k, Q);

            #pragma omp for
            for (size_t i = 1;                 i < m.rows()    - 1; i = i + 1)
            for (size_t j = 1 + ((i - 1) % 2); j < m.columns() - 1; j = j + 2)
                seidel_step(i, j, next, next, lamb, k, Q);

            #pragma omp single 
            {
                m.swap(next);
                iteration++;

                if ( (iteration % 100) == 0 ) {
                    helmholtz::write_matrix(next, std::to_string(iteration) + ".txt");
                }                
            }

            // std::cout << iteration << ": " << diff << std::endl;
            // if ( (iteration % 500) == 0 ) {
            //     helmholtz::writeMatrix(next, std::to_string(iteration) + ".txt");
            // }
        } while ( /*(diff > err) &&*/ (iteration < max_iterations) );
    }

    helmholtz::write_matrix(m, "out_conv.txt");
}
