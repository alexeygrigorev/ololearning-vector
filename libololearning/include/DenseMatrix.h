#include <cstddef>


class DenseVector {
public:
    DenseVector(size_t size);
    DenseVector(double* data, size_t size);
    ~DenseVector();

    double get(size_t i);
    void set(size_t i, double val);
    void swap(size_t i, size_t j);

    double norm2();

    double distance2(DenseVector other);
    double dot(DenseVector other);

    DenseVector subtract(DenseVector other, bool inplace);

    void printVector();

    size_t size();
    double* getData();

    DenseVector copy();
private:
    size_t _size;
    double* _data;
};


class DenseMatrix {
public:
    struct LUDecomposition {
        DenseMatrix *L;
        DenseMatrix *U;
    };


    DenseMatrix(size_t nrow, size_t ncol);
    DenseMatrix(double* data, size_t nrow, size_t ncol);
    ~DenseMatrix();

    static DenseMatrix eye(size_t n);
    static DenseMatrix ones(size_t nrow, size_t ncol);
    
    double get(size_t row, size_t col);
    void set(size_t row, size_t col, double val);

    double norm2();

    DenseVector getColumn(size_t col);
    DenseVector getRow(size_t row);
    void swapRows(size_t i, size_t j);

    double distance2(DenseMatrix other);
    DenseMatrix subtract(DenseMatrix other, bool inplace);

    DenseMatrix transpose(); 

    DenseVector vmult(DenseVector other);
    DenseMatrix mmult(DenseMatrix other);

    LUDecomposition lu();
    DenseVector solve(DenseVector b);
    DenseMatrix inverse();

    DenseMatrix copy();

    void printMatrix();

private:
    void init(double* data, size_t nrow, size_t ncol);
    size_t _nrow;
    size_t _ncol;
    double* _data;
};