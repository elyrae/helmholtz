#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <fstream>
#include <string>

#include "helmholtz.h"
#include "linear_matrix.h"

double  kappa(const double, const double) { return 1.0; }
double lambda(const double x, const double y) { return 1.0 / kappa(x, y); }
double source(const double x, const double y) {
    return   pow(M_E,-50*(pow(x,2) + pow(-0.5 + y,2)));
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
    // const double i = 1.0;

    // lg L = k(T)*log10(i*rho) + b(T)
    // return (k(T)*log10(D) + b(T))*log(10.0); // pow(k(T)*log10(i*rho) + b(T), 10.0);
    return exp( (k(T)*log10(D) + b(T))*log(10.0) );
}

void fill_data(LinearMatrix &kappa_equ, LinearMatrix &lambda_equ, LinearMatrix &Q_equ) 
{
    const std::string datafile = "dataLukin.dat";
    std::ifstream in(datafile);
    std::ofstream out("calculated.txt");

    const double sigma = 5.67E-5;             // Erg / s / cm^s / K^4
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

            if ((rho < 1.0E-16) || (T < 1000.0)) {
                // printf("(%d, %d): rho and T are too small. k = 0\n", i, j);
                kappa_equ(i, j) = 0.0;
                lambda_equ(i, j) = 1.0E12;
                Q_equ(i, j) = 0.0;

                out << x << " " << y << " " << 0.0 << "\n";
                continue;
            }
        
            kap = L_scale * cooling_function(rho, T) / (4.0*sigma*T*T*T*T);
            kappa_equ(i, j)  = 3.0 * kap;
            lambda_equ(i, j) = 1.0 / kap;
            Q_equ(i, j) = (3.0 / c) * L_scale * L_scale * cooling_function(rho, T); 

            // printf("(%d, %d): k = %.6f, q = %.6f\n", i, j, kap, Q_equ(i, j));
            out << x << " " << y << " " << cooling_function(rho, T) << "\n";
        }
    }
}

int main() {
    const size_t NX = 70;
    const size_t NY = 70;
    const double hx = 1.0 / double(NX - 1);
    const double hy = 1.0 / double(NY - 1);
    LinearMatrix m(NX, NY), k(NX, NY), lamb(NX, NY), Q(NX, NY);

    fill_data(k, lamb, Q);

    std::cout << cooling_function( exp(-3.707162E01), exp(6.907755E+00) ) << "\n";

    // double x = 0.0, y = 0.0;
    // for (size_t i = 0; i < m.rows();    ++i)
    // for (size_t j = 0; j < m.columns(); ++j) {
    //     x = i*hx;
    //     y = j*hy;

    //     k(i, j)    = kappa(x, y);
    //     lamb(i, j) = lambda(x, y);
    //     Q(i, j)    = source(x, y);
    // }

    // for (size_t i = 0; i < m.rows(); ++i) {
    //     // x = i*h;
    //     // y = i*h;

    //     m(0,            i) = 0.0;
    //     m(m.size() - 1, i) = 0.0;
    //     m(i,            0) = 0.0;
    //     m(i, m.size() - 1) = 0.0; 
    // }

    double start = omp_get_wtime();
    helmholtz::seidel_third_boundary(m, lamb, k, Q, 1.0E-4, 10000);
    double end = omp_get_wtime();

    printf("Time: %f s\n", end - start);
    return 0;
}
