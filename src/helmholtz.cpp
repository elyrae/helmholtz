#include "helmholtz.h"

#include <iostream>
#include <fstream>
//#include <iomanip>
#include <cmath>

void helmholtz::writeMatrix(const LinearMatrix& grid, const std::string& filename)
{
    std::ofstream out(filename);

    for (size_t i = 0; i < grid.rows(); ++i) {
        for (size_t j = 0; j < grid.columns(); ++j)
            out << grid(i, j) << " ";
        out << "\n"; 
    }
}

void helmholtz::jacobi(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
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
    double diff = 0.0, x = 0.0, y = 0.0, a_x_forw = 0.0, a_x_back = 0.0, a_y_forw = 0.0, a_y_back = 0.0;
    do {
        for (size_t i = 1; i < m.rows()    - 1; ++i)
        for (size_t j = 1; j < m.columns() - 1; ++j) {
            x = i*h;
            y = j*h;

            a_x_forw = lamb(x + h/2.0, y);
            a_x_back = lamb(x - h/2.0, y);
            a_y_forw = lamb(x, y + h/2.0);
            a_y_back = lamb(x, y - h/2.0);
            next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw 
                                      + m(i - 1, j)*a_x_back 
                                      + m(i, j + 1)*a_y_forw 
                                      + m(i, j - 1)*a_y_back) 
                         / (a_x_forw + a_x_back + a_y_forw + a_y_back + k(x, y)*h*h);
        }

        diff = 0.0;
        for (size_t i = 1; i < m.rows()    - 1; ++i)
        for (size_t j = 1; j < m.columns() - 1; ++j) {
            // if (fabs(next(i, j) - m(i, j)) > diff)
            //    diff = fabs(next(i, j) - m(i, j));
            m(i, j) = next(i, j);
        }

        // std::cout << iteration << ": " << diff << std::endl;
        // if ( ((iteration + 1) % 10) == 0 ) {
        //     helmholtz::writeMatrix(next, std::to_string(iteration + 1) + ".txt");
        // }

        iteration++;
        if ( (iteration % 100) == 0 ) {
            helmholtz::writeMatrix(next, std::to_string(iteration) + ".txt");
        }        
    } while ( /*(diff > err) &&*/ (iteration < max_iterations) );
}

double exact(const double x, const double y)
{
    return exp(-50.0 * ( (x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5) ));
}

void helmholtz::jacobiThirdBoundaryKind(LinearMatrix& m, 
                                        const ddFunction lamb, const ddFunction k, const ddFunction Q, 
                                        const double err, const size_t max_iterations) 
{
    LinearMatrix next(m.rows(), m.columns());

    const double h = 1.0 / double(m.size() - 1);

    size_t iteration = 0;
    double diff = 0.0, x = 0.0, y = 0.0, a_x_forw = 0.0, a_x_back = 0.0, a_y_forw = 0.0, a_y_back = 0.0;
    do {
        //=============================================
        // правая граница
        // for (size_t j = 1; j < m.columns() - 1; ++j) {
        //     double i = m.rows() - 1;
        //     x = i*h;
        //     y = j*h;

        //     a_x_forw = lamb(x + h/2.0, y        );
        //     a_x_back = lamb(x - h/2.0, y        );
        //     a_y_forw = lamb(x,         y + h/2.0);
        //     a_y_back = lamb(x,         y - h/2.0);

        //     next(i, j) = (Q(x, y)*h*h + m(i - 1, j)*(a_x_back + a_x_forw) + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
        //                / (a_x_forw*(1.0 + 3.0*h*k(x, y)) + a_y_forw + a_x_back + a_y_back + k(x, y)*h*h);
        // }

        //=============================================
        // левая граница
        // for (size_t j = 1; j < m.columns() - 1; ++j) {
        //     double i = 0.0;
        //     x = i*h;
        //     y = j*h;

        //     a_x_forw = lamb(x + h/2.0, y        );
        //     a_x_back = lamb(x - h/2.0, y        );
        //     a_y_forw = lamb(x,         y + h/2.0);
        //     a_y_back = lamb(x,         y - h/2.0);

        //     // next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back 
        //     //                           + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
        //     //            / (a_x_forw + a_y_forw + a_x_back + a_y_back + k(x, y)*h*h);

        //     next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*(a_x_forw + a_x_back) + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
        //                / (a_x_forw + a_y_forw + a_x_back*(1.0 + 3.0*h*k(x, y)) + a_y_back + k(x, y)*h*h);
        // }

        //=============================================
        // верхняя граница
        // for (size_t i = 1; i < m.rows() - 1; ++i) {
        //     double j = m.columns() - 1;
        //     x = i*h;
        //     y = j*h;

        //     a_x_forw = lamb(x + h/2.0, y        );
        //     a_x_back = lamb(x - h/2.0, y        );
        //     a_y_forw = lamb(x,         y + h/2.0);
        //     a_y_back = lamb(x,         y - h/2.0);

        //     next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back + m(i, j - 1)*(a_y_back + a_x_forw)) 
        //                / (a_x_forw + a_y_forw*(1.0 + 3.0*h*k(x, y)) + a_x_back + a_y_back + k(x, y)*h*h);

        // }

        //=============================================
        // нижняя граница
        // for (size_t i = 0; i < m.rows() - 1; ++i) {
        //     double j = 0;
        //     x = i*h;
        //     y = j*h;

        //     a_x_forw = lamb(x + h/2.0, y        );
        //     a_x_back = lamb(x - h/2.0, y        );
        //     a_y_forw = lamb(x,         y + h/2.0);
        //     a_y_back = lamb(x,         y - h/2.0);

        //     next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back + m(i, j + 1)*(a_y_forw + a_y_back))
        //                / (a_x_forw + a_y_forw + a_x_back + a_y_back*(1.0 + 3.0*h*k(x, y)) + k(x, y)*h*h);
        // }

        // helmholtz::writeMatrix(next, "right_left_up_down.txt");
        
        // пересчет
        // x = 0.0;
        // y = 0.0;
        
        // a_x_forw = lamb(x + h/2.0, y        );
        // a_x_back = lamb(x - h/2.0, y        );
        // a_y_forw = lamb(x,         y + h/2.0);
        // a_y_back = lamb(x,         y - h/2.0);
        // next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back 
        //                           + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
        //            / (a_x_forw + a_y_forw + a_x_back + a_y_back + k(x, y)*h*h);
        
        // пересчет
        // x = i*h;
        // y = j*h;
        
        // a_x_forw = lamb(x + h/2.0, y        );
        // a_x_back = lamb(x - h/2.0, y        );
        // a_y_forw = lamb(x,         y + h/2.0);
        // a_y_back = lamb(x,         y - h/2.0);
        // next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back 
        //                           + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
        //            / (a_x_forw + a_y_forw + a_x_back + a_y_back + k(x, y)*h*h);
        
        // пересчет
        // x = i*h;
        // y = j*h;
        
        // a_x_forw = lamb(x + h/2.0, y        );
        // a_x_back = lamb(x - h/2.0, y        );
        // a_y_forw = lamb(x,         y + h/2.0);
        // a_y_back = lamb(x,         y - h/2.0);
        // next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back 
        //                           + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
        //            / (a_x_forw + a_y_forw + a_x_back + a_y_back + k(x, y)*h*h);
        
        // пересчет
        // x = i*h;
        // y = j*h;
        
        // a_x_forw = lamb(x + h/2.0, y        );
        // a_x_back = lamb(x - h/2.0, y        );
        // a_y_forw = lamb(x,         y + h/2.0);
        // a_y_back = lamb(x,         y - h/2.0);
        // next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back 
        //                           + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
        //            / (a_x_forw + a_y_forw + a_x_back + a_y_back + k(x, y)*h*h);                                      

        //=============================================
        // пересчет внутренних точек сетки
        // for (size_t i = 1; i < m.rows()    - 1; ++i)
        // for (size_t j = 1; j < m.columns() - 1; ++j) {
        //     x = i*h;
        //     y = j*h;

        //     a_x_forw = lamb(x + h/2.0, y        );
        //     a_x_back = lamb(x - h/2.0, y        );
        //     a_y_forw = lamb(x,         y + h/2.0);
        //     a_y_back = lamb(x,         y - h/2.0);

        //     next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back 
        //                               + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
        //                / (a_x_forw + a_x_back + a_y_forw + a_y_back + k(x, y)*h*h);
        // }

        for (size_t i = 0; i < m.rows()   ; ++i)
        for (size_t j = 0; j < m.columns(); ++j) {
            x = i*h;
            y = j*h;

            a_x_forw = lamb(x + h/2.0, y        );
            a_x_back = lamb(x - h/2.0, y        );
            a_y_forw = lamb(x,         y + h/2.0);
            a_y_back = lamb(x,         y - h/2.0);

            next(i, j) = (  Q(x, y)*h*h 
                          + m(i + 1, j)*(a_x_forw + a_x_back*double(i == 0                ))*double(i != (m.rows() - 1)   ) 
                          + m(i - 1, j)*(a_x_back + a_x_forw*double(i == (m.rows() - 1)   ))*double(i != 0                ) 
                          + m(i, j + 1)*(a_y_forw + a_y_back*double(j == 0                ))*double(j != (m.columns() - 1)) 
                          + m(i, j - 1)*(a_y_back + a_y_forw*double(j == (m.columns() - 1)))*double(j != 0                )) 
                         / 
                         (  a_x_forw*(1.0 + 3.0*h*k(x, y)*double(i == (m.rows() - 1))   ) 
                          + a_x_back*(1.0 + 3.0*h*k(x, y)*double(i == 0)                ) 
                          + a_y_forw*(1.0 + 3.0*h*k(x, y)*double(j == (m.columns() - 1))) 
                          + a_y_back*(1.0 + 3.0*h*k(x, y)*double(j == 0)                ) 
                          + k(x, y)*h*h);
        }

        // diff = 0.0;
        for (size_t i = 1; i < m.rows()    - 1; ++i)
        for (size_t j = 1; j < m.columns() - 1; ++j) {
            // if (fabs(next(i*h, j*h) - exact(x, y)) > diff)
            //   diff = fabs(next(i*h, j*h) - exact(x, y));

            m(i, j) = next(i, j);
        }

        iteration++;
        // std::cout << iteration << ": " << diff << std::endl;
        if ( (iteration % 500) == 0 ) {
            helmholtz::writeMatrix(next, std::to_string(iteration) + ".txt");
        }

        // printf("%.15f\n", diff);

        //iteration++;
    } while ( /*(diff > err) &&*/ (iteration < max_iterations) );
}

// void helmholtz::seidel(LinearMatrix& m, const ddFunction lamb, const ddFunction k, const ddFunction Q, 
//                        const double err, const size_t max_iterations) 
// {
//     LinearMatrix next(m.rows(), m.columns());
//     for (size_t i = 0; i < m.rows(); ++i) {
//         next(0,            i) = m(0,            i);
//         next(m.size() - 1, i) = m(m.size() - 1, i);
//         next(i,            0) = m(i,            0);
//         next(i, m.size() - 1) = m(i, m.size() - 1); 
//     }

//     const double h = 1.0 / double(m.size() - 1);

//     size_t iteration = 0;
//     double diff = 0.0, x = 0.0, y = 0.0, a_x_forw = 0.0, a_x_back = 0.0, a_y_forw = 0.0, a_y_back = 0.0;    
//     do {
//         for (size_t i = 1;           i < m.rows()    - 1; i = i + 1)
//         for (size_t j = 1 + (i % 2); j < m.columns() - 1; j = j + 2) {
//             x = i*h;
//             y = j*h;

//             a_x_forw = lamb(x + h/2.0, y);
//             a_x_back = lamb(x - h/2.0, y);
//             a_y_forw = lamb(x, y + h/2.0);
//             a_y_back = lamb(x, y - h/2.0);

//             next(i, j) = (Q(x, y)*h*h + m(i + 1, j)*a_x_forw + m(i - 1, j)*a_x_back 
//                                       + m(i, j + 1)*a_y_forw + m(i, j - 1)*a_y_back) 
//                        / (a_x_forw + a_x_back + a_y_forw + a_y_back + k(x, y)*h*h);   
//         }

//         for (size_t i = 1;                 i < m.rows()    - 1; i = i + 1)
//         for (size_t j = 1 + ((i - 1) % 2); j < m.columns() - 1; j = j + 2) {
//             x = i*h;
//             y = j*h;

//             a_x_forw = lamb(x + h/2.0, y);
//             a_x_back = lamb(x - h/2.0, y);
//             a_y_forw = lamb(x, y + h/2.0);
//             a_y_back = lamb(x, y - h/2.0);

//             next(i, j) = (Q(x, y)*h*h + next(i + 1, j)*a_x_forw + next(i - 1, j)*a_x_back 
//                                       + next(i, j + 1)*a_y_forw + next(i, j - 1)*a_y_back) 
//                        / (a_x_forw + a_x_back + a_y_forw + a_y_back + k(x, y)*h*h);
//         }

//         diff = 0.0;
//         for (size_t i = 1; i < m.rows()    - 1; ++i)
//         for (size_t j = 1; j < m.columns() - 1; ++j) {
//             if (fabs(next(i, j) - m(i, j)) > diff)
//                 diff = fabs(m(i, j) - next(i, j));
//             m(i, j) = next(i, j);
//         }

//         std::cout << iteration << ": " << diff << std::endl;
//         if ( ((iteration + 1) % 10) == 0 ) {
//             helmholtz::writeMatrix(next, std::to_string(iteration + 1) + ".txt");
//         }

//         iteration++;
//     } while ( (diff > err) && (iteration < max_iterations) );
// }
