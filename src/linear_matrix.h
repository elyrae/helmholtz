#include <vector>

class LinearMatrix {
    size_t matrix_rows;    
    size_t matrix_columns;
    std::vector<double> v;
public:
    LinearMatrix(const size_t rows, const size_t columns) : matrix_rows(rows), matrix_columns(columns), v(rows*columns, 0.0) {}

          double& operator()(const size_t i, const size_t j)       { return v[i*matrix_columns + j]; }
    const double& operator()(const size_t i, const size_t j) const { return v[i*matrix_columns + j]; }
    size_t rows()    const { return matrix_rows; }
    size_t size()    const { return rows(); }
    size_t columns() const { return matrix_columns; }

    void fillZero() {
        for (size_t i = 0; i < rows()*columns(); ++i) 
            v[i] = 0.0;   
    }
};
