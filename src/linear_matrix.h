#ifndef LINEAR_MATRIX_H
#define LINEAR_MATRIX_H

#include <vector>
#include <cstddef>

class LinearMatrix {
    size_t matrix_rows;    
    size_t matrix_columns;
    std::vector<double> v;
public:
    LinearMatrix(const size_t rows, const size_t columns) : matrix_rows(rows), matrix_columns(columns), v(rows*columns) {}

    LinearMatrix &operator=(const LinearMatrix &that) 
    {
        matrix_rows    = that.rows();
        matrix_columns = that.columns();
        v = that.vect();

        return *this;
    }

    ~LinearMatrix() {}

          double& operator()(const size_t i, const size_t j)       { return  v[i*matrix_columns + j]; }
    const double& operator()(const size_t i, const size_t j) const { return  v[i*matrix_columns + j]; }
          double* operator[](const size_t i)                       { return &v[i*matrix_columns];     }
    const double* operator[](const size_t i)                 const { return &v[i*matrix_columns];     }    
    size_t rows()    const { return matrix_rows;    }
    size_t size()    const { return rows();         }
    size_t columns() const { return matrix_columns; }

    double *data() { return v.data(); }

          std::vector<double>& vect()       { return v; }
    const std::vector<double>& vect() const { return v; }

    void resize(const size_t rows, const size_t columns) 
    {
        matrix_rows = rows;
        matrix_columns = columns;
        v.resize(rows*columns);
    }

    void clear() { std::fill(v.begin(), v.end(), 0); }
    
    void swap(LinearMatrix& second) { v.swap(second.vect()); }
};

#endif