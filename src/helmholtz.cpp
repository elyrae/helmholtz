#include "helmholtz.h"

#include <iostream>
#include <fstream>
#include <cmath>

const int THREADS = 4;

void helmholtz::write_matrix(const LinearMatrix &grid, const std::string &filename)
{
    std::ofstream out(filename);
    for (size_t i = 0; i < grid.rows(); ++i) {
        for (size_t j = 0; j < grid.columns(); ++j)
            out << grid(i, j) << " ";
        out << "\n"; 
    }
}

void helmholtz::jacobi_first_boundary(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                                      const double err, const size_t max_iterations) 
{
    LinearMatrix next = m; 
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

        m.swap(next);
        iteration++;
        if ( (iteration % 100) == 0 ) {
            helmholtz::write_matrix(next, std::to_string(iteration) + ".txt");
        }        
    } while (iteration < max_iterations);
}

void seidel_step(const size_t i, const size_t j, LinearMatrix &next, const LinearMatrix &m, 
                 const LinearMatrix &lamb, const LinearMatrix &k, const LinearMatrix &Q)
{
    const double hx = 1.0 / double(m.rows()    - 1);
    const double hy = 1.0 / double(m.columns() - 1);
    const double hy_hx = hy / hx;
    const double hx_hy = hx / hy;

    const double    k_ij = k(i, j);
    const double lamb_ij = lamb(i, j);
    const double    q_ij = Q(i, j);

    const double a_x_forw = (lamb(i + 1, j) + lamb(i, j)) / 2.0;
    const double a_x_back = (lamb(i - 1, j) + lamb(i, j)) / 2.0;
    const double a_y_forw = (lamb(i, j + 1) + lamb(i, j)) / 2.0;
    const double a_y_back = (lamb(i, j - 1) + lamb(i, j)) / 2.0;

    // =========================================================== 
    // внутренняя часть области
    if ((0 < i) && (i < (m.rows() - 1)) && (0 < j) && (j < (m.columns() - 1)))
        next(i, j) = (q_ij*hy*hy + m(i + 1, j)*a_x_forw*hy_hx*hy_hx 
                                 + m(i - 1, j)*a_x_back*hy_hx*hy_hx 
                                 + m(i, j + 1)*a_y_forw 
                                 + m(i, j - 1)*a_y_back) 
                      / ((a_x_forw + a_x_back)*hy_hx*hy_hx + a_y_forw + a_y_back + k_ij*hy*hy);

    // ===========================================================
    // левая граница без угловых точек
    if ((i == 0) && (0 < j) && (j < (m.columns() - 1)))
        next(i, j) = (q_ij*hy*hy*0.5 + m(1, j    )*a_x_forw*hy_hx*hy_hx 
                                     + m(0, j + 1)*a_y_forw*0.5 
                                     + m(0, j - 1)*a_y_back*0.5) 
                     / (a_x_forw*hy_hx*hy_hx + 0.5*k_ij*lamb_ij*hy*hy_hx + 0.5*a_y_forw + 0.5*a_y_back + k_ij*hy*hy*0.5);

    // правая граница без угловых точек
    if ((i == (m.rows() - 1)) && (0 < j) && (j < (m.columns() - 1)))
        next(i, j) = (q_ij*hy*hy*0.5 + m(i - 1, j)*a_x_back*hy_hx*hy_hx 
                                     + m(i, j + 1)*a_y_forw*0.5 
                                     + m(i, j - 1)*a_y_back*0.5) 
                     / (a_x_back*hy_hx*hy_hx + 0.5*k_ij*lamb_ij*hy*hy_hx + 0.5*a_y_forw + 0.5*a_y_back + k_ij*hy*hy*0.5);

    // верхняя граница без угловых точек
    if ((0 < i) && (i < (m.rows() - 1)) && (j == (m.columns() - 1)))
        next(i, j) = (q_ij*hx*hx*0.5 + m(i + 1, j)*a_x_forw*0.5 
                                     + m(i - 1, j)*a_x_back*0.5 
                                     + m(i, j - 1)*a_y_back*hx_hy*hx_hy) 
                     / (0.5*a_x_forw + 0.5*a_x_back + 0.5*k_ij*lamb_ij*hx*hx_hy + a_y_back*hx_hy*hx_hy + k_ij*hx*hx*0.5);

    // нижняя граница без угловых точек
    if ((0 < i) && (i < (m.rows() - 1)) && (j == 0))
        next(i, j) = (q_ij*hx*hx*0.5 + m(i + 1, 0)*a_x_forw*0.5 
                                     + m(i - 1, 0)*a_x_back*0.5 
                                     + m(i,     1)*a_y_forw*hx_hy*hx_hy) 
                     / (0.5*a_x_forw + 0.5*a_x_back + a_y_forw*hx_hy*hx_hy + 0.5*k_ij*lamb_ij*hx*hx_hy + k_ij*hx*hx*0.5);

    // ===========================================================                  
    // левый верхний угол +
    if ((i == 0) && (j == (m.columns() - 1)))
        next(i, j) = (q_ij*hy*hy*0.5 + m(1, j    )*a_x_forw*hy_hx*hy_hx 
                                     + m(0, j - 1)*a_y_back)
                     / (a_x_forw*hy_hx*hy_hx + 0.5*k_ij*lamb_ij*hy + 0.5*k_ij*lamb_ij*hy*hy_hx + a_y_back + k_ij*hy*hy*0.5);

    // левый нижний угол +
    if ((i == 0) && (j == 0))
        next(i, j) = (q_ij*hy*hy*0.5 + m(1, 0)*a_x_forw*hy_hx*hy_hx  
                                     + m(0, 1)*a_y_forw) 
                     / (a_x_forw*hy_hx*hy_hx + 0.5*k_ij*lamb_ij*hy + a_y_forw + 0.5*k_ij*lamb_ij*hy*hy_hx + k_ij*hy*hy*0.5);

    // правый верхний угол +++
    if ((i == (m.rows() - 1)) && (j == (m.columns() - 1)))
        next(i, j) = (q_ij*hy*hy*0.5 + m(i - 1, j)*a_x_back*hy_hx*hy_hx 
                                     + m(i, j - 1)*a_y_back) 
                     / (0.5*k_ij*lamb_ij*hy + a_x_back*hy_hx*hy_hx + 0.5*k_ij*lamb_ij*hy*hy_hx + a_y_back + k_ij*hy*hy*0.5); 

    // правый нижний угол
    if ((i == (m.rows() - 1)) && (j == 0)) 
        next(i, j) = (q_ij*hy*hy*0.5 + m(i - 1, 0)*a_x_back*hy_hx*hy_hx 
                                     + m(i,     1)*a_y_forw) 
                     / (0.5*k_ij*lamb_ij*hy + a_x_back*hy_hx*hy_hx + a_y_forw + 0.5*k_ij*lamb_ij*hy*hy_hx + k_ij*hy*hy*0.5);
}

void helmholtz::seidel_third_boundary(LinearMatrix& m, 
                                      const LinearMatrix &lamb, const LinearMatrix &k, const LinearMatrix &Q, 
                                      const double err, const size_t max_iterations) 
{
    LinearMatrix next(m.rows(), m.columns());

    size_t iteration = 0;
    do {
        // #pragma omp for
        for (size_t i = 0;           i < m.rows()   ; i = i + 1)
        for (size_t j = 0 + (i % 2); j < m.columns(); j = j + 2)
            seidel_step(i, j, next, m, lamb, k, Q);
    
        // #pragma omp for
        for (size_t i = 0;                 i < m.rows()   ; i = i + 1)
        for (size_t j = 0 + ((i + 1) % 2); j < m.columns(); j = j + 2)
            seidel_step(i, j, next, next, lamb, k, Q);
    
        m.swap(next);
        iteration++;
    } while (iteration < max_iterations);
}
