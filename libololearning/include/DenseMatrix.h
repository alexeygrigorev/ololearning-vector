#include <cstddef>


class DenseVector {
public:
    DenseVector(size_t size);
    DenseVector(double* data, size_t size);
    double get(size_t i);
    void set(size_t i, double val);
    double norm2();

    double distance2(DenseVector other);
    double dot(DenseVector other);

    DenseVector subtract(DenseVector other, bool inplace);

    void printVector();
private:
    size_t size;
    double* data;
};


class DenseMatrix {
public: 
    DenseMatrix(size_t nrow, size_t ncol);
    static DenseMatrix eye(size_t n);
    static DenseMatrix ones(size_t nrow, size_t ncol);

    DenseMatrix(double* data, size_t nrow, size_t ncol);
    double get(size_t row, size_t col);
    void set(size_t row, size_t col, double val);

    double norm2();

    DenseVector getColumn(size_t col);
    DenseVector getRow(size_t row);

    double distance2(DenseMatrix other);
    DenseMatrix subtract(DenseMatrix other, bool inplace);

    DenseMatrix transpose(); 

    DenseVector vmult(DenseVector other);
    DenseMatrix mmult(DenseMatrix other);

    void printMatrix();

private:
    void init(double* data, size_t nrow, size_t ncol);
    size_t nrow;
    size_t ncol;
    double* data;
};