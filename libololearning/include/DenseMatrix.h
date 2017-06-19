#include <cstddef>

#pragma once
#include "DenseVector.h"

class DenseMatrix;

struct LUDecomposition {
    DenseMatrix *P;
    DenseMatrix *L;
    DenseMatrix *U;
};


struct QRDecomposition {
    DenseMatrix *QT;
    DenseMatrix *R;
};


class DenseMatrix {
public:

    DenseMatrix(size_t nrow, size_t ncol);
    DenseMatrix(float* data, size_t nrow, size_t ncol);
    ~DenseMatrix();

    static DenseMatrix* eye(size_t n);
    static DenseMatrix* ones(size_t nrow, size_t ncol);
    
    float get(size_t row, size_t col);
    void set(size_t row, size_t col, float val);

    float norm2();

    DenseVector* getColumn(size_t col);
    DenseVector* getRow(size_t row);
    float* getRowData(size_t row);

    void swapRows(size_t i, size_t j);

    float distance2(DenseMatrix *other);
    DenseMatrix* subtract(DenseMatrix *other, bool inplace);

    DenseMatrix* transpose(); 

    DenseVector* vmult(DenseVector *other);
    DenseMatrix* mmult(DenseMatrix *other);

    LUDecomposition lu();

    DenseVector* solve(DenseVector *b);
    DenseMatrix* solveMatrix(DenseMatrix *B);

    DenseMatrix* inverseGaussJordan();
    DenseMatrix* inverseLU();
    DenseMatrix* inverseQR();
    DenseMatrix* inverse();

    DenseMatrix* orthonormalize();

    QRDecomposition qr();

    DenseMatrix* copy();

    size_t numRows();
    size_t numCols();

    void printMatrix();

private:
    DenseMatrix(float* data, size_t nrow, size_t ncol, bool external);
    void init(float* data, size_t nrow, size_t ncol, bool external);

    size_t _nrow;
    size_t _ncol;
    float* _data;
    bool _external;
};