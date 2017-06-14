#include <cstddef>


class DenseVector {
public:
    DenseVector(size_t size);
    DenseVector(float* data, size_t size);
    ~DenseVector();

    float get(size_t i);
    void set(size_t i, float val);
    void swap(size_t i, size_t j);

    float norm2();

    float distance2(DenseVector other);
    float dot(DenseVector other);

    DenseVector subtract(DenseVector other, bool inplace);

    void printVector();

    size_t size();
    float* getData();

    DenseVector copy();
private:
    size_t _size;
    float* _data;
};

class DenseMatrix;

struct LUDecomposition {
    DenseMatrix *P;
    DenseMatrix *L;
    DenseMatrix *U;
};


class DenseMatrix {
public:

    DenseMatrix(size_t nrow, size_t ncol);
    DenseMatrix(float* data, size_t nrow, size_t ncol);
    ~DenseMatrix();

    static DenseMatrix eye(size_t n);
    static DenseMatrix ones(size_t nrow, size_t ncol);
    
    float get(size_t row, size_t col);
    void set(size_t row, size_t col, float val);

    float norm2();

    DenseVector getColumn(size_t col);
    DenseVector getRow(size_t row);
    void swapRows(size_t i, size_t j);

    float distance2(DenseMatrix other);
    DenseMatrix subtract(DenseMatrix other, bool inplace);

    DenseMatrix transpose(); 

    DenseVector vmult(DenseVector other);
    DenseMatrix mmult(DenseMatrix other);

    LUDecomposition lu();

    DenseVector gaussJordanEliminationVector(DenseVector b);
    DenseMatrix gaussJordanEliminationMatrix(DenseMatrix B);

    DenseVector solve(DenseVector b);
    DenseMatrix solveMatrix(DenseMatrix B);
    DenseMatrix inverse();

    DenseMatrix copy();

    size_t numRows();
    size_t numCols();

    void printMatrix();

private:
    void init(float* data, size_t nrow, size_t ncol);
    size_t _nrow;
    size_t _ncol;
    float* _data;
};