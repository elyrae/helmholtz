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
    return   pow(M_E,-50*(pow(x,2) + pow(y,2))) - pow(M_E,-50*(pow(-1 + x,2) + pow(-1 + y,2)));
}

double k(const double T)
{
    if (log10(T) < 3.65)
        return 1.0;
    else if (log10(T) < 4.15)
        return 2.0*log10(T) - 6.3;
    else
        return 2.0;
}

double b(const double T)
{
    if (log10(T) < 3.325)
        return 86.154*log10(T) - 286.46;
    else if (log10(T) < 4.1)
        return 33.226*log10(T) - 110.48;
    else if (log10(T) < 5.25)
        return -1.957*log10(T) + 33.77;
    else
        return 0.667*log10(T) + 20.0;
}

double cooling_function(const double rho, const double T)
{
    const double i = 1.0;

    // lg L = k(T)*log10(i*rho) + b(T)
    return pow(k(T)*log10(i*rho) + b(T), 10.0);
}

void fill_data(LinearMatrix &kappa, LinearMatrix &Q) 
{
    const std::string datafile = "dataLukin.dat";
    std::ifstream in(datafile);

    const double sigma   = 5.67E-5;  // Erg / s / cm^s / K^4
    const double R_sun   = 6.957E10; // cm
    const double c       = 299792458.0 * 100.0; // cm / s 

    const double L_scale = R_sun * 0.557727 * 1.38;

    double x = 0.0, y = 0.0, rho = 0.0, T = 0.0, energy = 0.0;
    for (size_t i = 0; i < 70; ++i)
    for (size_t j = 0; j < 70; ++j) {
        in >> x;
        in >> y;
        in >> rho;
        in >> T;
        in >> energy;

        if (rho < 1.0E-16) {
            printf("(%d, %d): rho is too small. kappa = 0\n", i, j);
            kappa(i, j) = 0.0;

            continue;
        }

        if (T < 1000.0) {
            printf("(%d, %d): T is too small. T = 0\n", i, j);
            kappa(i, j) = 0.0;

            continue;
        }
        
        kappa(i, j) = L_scale * cooling_function(rho, T) / (4.0*sigma*T*T*T*T);
        printf("(%d, %d): k = %.12f\n", i, j, kappa(i, j));
    }
}

int main() {
    const size_t NX = 70;
    const size_t NY = 70;
    const double hx = 1.0 / double(NX - 1);
    const double hy = 1.0 / double(NY - 1);
    LinearMatrix m(NX, NY), k(NX, NY), lamb(NX, NY), Q(NX, NY);

    fill_data(k, Q);

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

    // double start = omp_get_wtime();
    // helmholtz::seidel_third_boundary(m, lamb, k, Q, 1.0E-4, 10000);
    // double end = omp_get_wtime();

    // printf("Time: %f s\n", end - start);
    return 0;
}
