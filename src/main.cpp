#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <fstream>
#include <string>

#include "helmholtz.h"
#include "linear_matrix.h"

double kp(const double x, const double y) {
    return 1.0E-6;
}

double source(const double x, const double y) {
    return (-2*pow(-1 + y,2)*pow(y,2) + 12*x*pow(-1 + y,2)*pow(y,2) + pow(x,3)*(4 - 2*(-4 + y)*(-1 + y)*y*(3 + y)) + 
   pow(x,4)*(-2 + (-4 + y)*(-1 + y)*y*(3 + y)) + pow(x,2)*(-2 - (-1 + y)*y*(12 + 11*(-1 + y)*y)) );
    // return (-4799 + 10000*x - 10000*pow(x,2) + 10000*y - 10000*pow(y,2))/pow(M_E,50*(pow(-0.5 + x,2) + pow(-0.5 + y,2))); 
    // pow(M_E,-50*(pow(x,2) + pow(-0.5 + y,2)));
}

double k(const double T)
{
    const double logT = log10(T);
    if (logT < 3.65)
        return 1.0;
    else if (logT < 4.15)
        return 2.0*logT - 6.3;
    else
        return 2.0;
}

double b(const double T)
{
    const double logT = log10(T);
    if (logT < 3.325)
        return 86.154*logT - 286.46;
    else if (logT < 4.1)
        return 33.226*logT - 110.48;
    else if (logT < 5.25)
        return -1.957*logT + 33.77;
    else
        return  0.667*logT + 20.0;
}

double cooling_function(const double D, const double T)
{
    return exp( (k(T)*log10(D) + b(T))*log(10.0) );
}

void fill_data(LinearMatrix &kappa_equ, LinearMatrix &lambda_equ, LinearMatrix &Q_equ) 
{
    const std::string datafile = "dataLukin.dat";
    std::ifstream in(datafile);
    
    std::ofstream out_kappa("kappa.txt");
    std::ofstream out_lambda("lambda.txt");
    std::ofstream out_Q("Q.txt");

    const double sigma = 5.67E-5;             // Erg / s / cm^2 / K^4
    const double R_sun = 6.957E10;            // cm
    const double c     = 299792458.0 * 100.0; // cm / s 

    const double L_scale = R_sun * 0.557727 * 1.38;

    double x = 0.0, y = 0.0, rho = 0.0, T = 0.0, energy = 0.0, kap = 0.0;
    for (size_t i = 0; i < 70; ++i) {
        for (size_t j = 0; j < 70; ++j) {
            in >> x;
            in >> y;
            in >> rho;
            in >> T;
            in >> energy;

            if ((rho < 1.0E-14) || (T < 1000.0)) {
                kap = 1.0E-15;

                kappa_equ(i, j) = 3.0 * kap;
                lambda_equ(i, j) = 1.0E15;
                Q_equ(i, j)     = 3.0 * kap * L_scale * energy; // 0.0;
            } 
            else {
                kap = L_scale * cooling_function(rho, T) / (4.0*sigma*T*T*T*T);

                kappa_equ(i, j)  =  3.0 * kap;
                lambda_equ(i, j) = (1.0 / kap);
                Q_equ(i, j) = 3.0 * kap * L_scale * energy; // (3.0 / c) * L_scale * L_scale * cooling_function(rho, T); 
            }    

            out_kappa  << x << " " << y << " " << kappa_equ(i, j)  << "\n";
            out_lambda << x << " " << y << " " << lambda_equ(i, j) << "\n";
            out_Q      << x << " " << y << " " << Q_equ(i, j) << "\n";
        }
    }
}

double sol(const double x, const double y)
{
    return x*x*y*y*(1.0 - x)*(1.0 - x)*(1.0 - y)*(1.0 - y); 
    // pow(M_E, -50*(pow(-0.5 + x,2) + pow(-0.5 + y,2)));
}

void show_error(const LinearMatrix &m, const size_t NX, const size_t NY, const double hx, const double hy)
{
    double max = 0.0, x = 0.0, y = 0.0;
    for (size_t i = 0; i < NX; ++i)
    for (size_t j = 0; j < NY; ++j) {
        x = i*hx;
        y = j*hy;

        if ( fabs(m(i, j) - sol(x, y)) > max )
            max = fabs(m(i, j) - sol(x, y));
    }

    std::cout << "err: " << max << "\n";
}

int main() {
    const size_t NX = 70;
    const size_t NY = 70;
    const double hx = 1.0 / double(NX - 1);
    const double hy = 1.0 / double(NY - 1);
    LinearMatrix m(NX, NY), k(NX, NY), lamb(NX, NY), Q(NX, NY);

    fill_data(k, lamb, Q);

    double kap = 0.0, x = 0.0, y = 0.0;
    for (size_t i = 0; i < NX; ++i)
    for (size_t j = 0; j < NY; ++j) {
        x = i*hx;
        y = j*hy;

        // kap = kp(x, y);
        k(i, j)    = 1.0; // 3.0 * kap;
        lamb(i, j) = 1.0; //1.0 / kap;
        Q(i, j) = source(x, y);

        // Q(NX/2, NY/2) = 3.0 * kp(NX/2, NY/2) * 1.0;
    }

    double start = omp_get_wtime();
    helmholtz::seidel_third_boundary(m, lamb, k, Q, 1.0E-6, 500000);
    double end = omp_get_wtime();

    // show_error(m, NX, NY, hx, hy);
    // const std::string filename = "post";
    // write_solution_vtk(filename + ".txt", m, NX, NY);

    printf("Time: %f s\n", end - start);
    return 0;
}
