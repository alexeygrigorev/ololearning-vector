#include <cstddef>
#include <cstring>
#include <iostream>
#include <cmath>

#include <stdexcept>
#include <cassert>

#include "DenseVector.h"
#include "DenseMatrix.h"
#include "math_utils.h"

using namespace std;
using namespace olomath;

void DenseMatrix::init(float* data, size_t nrow, size_t ncol, bool external) {
    this->_nrow = nrow;
    this->_ncol = ncol;
    this->_data = data;
    this->_external = external;
}

DenseMatrix::DenseMatrix(float* data, size_t nrow, size_t ncol) {
    this->init(data, nrow, ncol, true);
}

DenseMatrix::DenseMatrix(float* data, size_t nrow, size_t ncol, bool external) {
    this->init(data, nrow, ncol, external);
}


DenseMatrix::DenseMatrix(size_t nrow, size_t ncol) {
    float* data = zeroMatrix(nrow, ncol);
    this->init(data, nrow, ncol, false);
}

DenseMatrix::~DenseMatrix() {
    cout << "Dense Matrix destructor" << endl;
    if (!this->_external) {
        cout << "removing the array" << endl;
        delete[] this->_data;
        this->_data = nullptr;
    }
}


DenseMatrix* DenseMatrix::eye(size_t n) {
    DenseMatrix *m = new DenseMatrix(n, n);

    for (size_t i = 0; i < n; i++) {
        m->set(i, i, 1);
    }

    return m;
}


float DenseMatrix::get(size_t row, size_t col) {
    assert(row < this->_nrow);
    assert(col < this->_ncol);
    size_t ncol = this->_ncol;
    return this->_data[row * ncol + col];
}

void DenseMatrix::set(size_t row, size_t col, float val) {
    assert(row < this->_nrow);
    assert(col < this->_ncol);
    size_t ncol = this->_ncol;
    float* data = this->_data;
    data[row * ncol + col] = val;
}

DenseVector* DenseMatrix::getColumn(size_t col) {
    assert(col < this->_ncol);
    size_t nrow = this->_nrow;

    float* colData = new float[nrow];

    for (size_t r = 0; r < nrow; r++) {
        colData[r] = this->get(r, col);
    }

    return new DenseVector(colData, nrow, false);
}

float* DenseMatrix::getRowData(size_t row) {
    assert(row < this->_nrow);
    float* data = this->_data;
    size_t ncol = this->_ncol;
    return &data[row * ncol];
}

DenseVector* DenseMatrix::getRow(size_t row) {
    size_t ncol = this->_ncol;
    float* rowData = getRowData(row);
    return new DenseVector(rowData, ncol, true);
}

void DenseMatrix::swapRows(size_t i, size_t j) {
    assert(i < this->_nrow);
    assert(j < this->_nrow);

    if (i == j) {
        return;
    }

    float* data = this->_data;
    size_t ncol = this->_ncol;

    float* rowI = &data[i * ncol];
    float* rowJ = &data[j * ncol];
    float* tmp = new float[ncol];

    size_t size = ncol * sizeof(float);
    memcpy(tmp, rowI, size);
    memcpy(rowI, rowJ, size);
    memcpy(rowJ, tmp, size);
    delete[] tmp;
}


DenseMatrix* DenseMatrix::subtract(DenseMatrix *other, bool inplace) {
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;
    assert(ncol == other->_ncol);
    assert(nrow == other->_nrow);

    size_t size = ncol * nrow;
    float* newData = copyOrSame(this->_data, size, inplace);
    float* data = this->_data;
    float* oData = other->_data;

    for (size_t i = 0; i < size; i++) {
        newData[i] = data[i] - oData[i];
    }

    if (inplace == true) {
        return this;
    } else {
        return new DenseMatrix(newData, nrow, ncol, false);
    }
}

float DenseMatrix::norm2() {
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;
    size_t size = ncol * nrow;

    float* data = this->_data;
    return arrayNorm2(data, size);
}

float DenseMatrix::distance2(DenseMatrix *other) {
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;
    assert(ncol == other->_ncol);
    assert(nrow == other->_nrow);

    size_t size = ncol * nrow;

    float* data = this->_data;
    float* oData = other->_data;

    float total = 0.0;
    for (size_t i = 0; i < size; i++) {
        float d = data[i] - oData[i];
        total = total + d * d;
    }

    return total;
}

DenseMatrix* DenseMatrix::transpose() {
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;

    DenseMatrix *t = new DenseMatrix(ncol, nrow);
    t->_external = false;

    for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            float val = this->get(i, j);
            t->set(j, i, val);
        }
    }

    return t;
}


DenseVector* DenseMatrix::vmult(DenseVector *vec) {
    size_t nrow = this->_nrow;
    size_t ncol = this->_ncol;
    assert(this->_ncol == vec->size());

    float* result = new float[nrow];

    float* vecData = vec->getData();

    for (size_t i = 0; i < nrow; i++) {
        float *row = this->getRowData(i);
        result[i] = arrayDot(row, vecData, ncol);
    }

    return new DenseVector(result, nrow);
}

DenseMatrix* DenseMatrix::mmult(DenseMatrix *other) {
    assert(this->_ncol == other->_nrow);
    size_t nrowRes = this->_nrow;
    size_t ncol = this->_ncol;
    size_t ncolRes = other->_ncol;

    size_t size = nrowRes * ncolRes;
    float *thisMatrix = this->_data;

    float *resBuffer = new float[nrowRes];
    size_t resBufferSize = nrowRes * sizeof(float);
    float *resultData = new float[size];
    float *dataPtr = resultData;

    DenseMatrix* o = other->transpose();

    for (size_t i = 0; i < ncolRes; i++) {
        float *oCol = o->getRowData(i);
        matrixVectorProduct(thisMatrix, oCol, resBuffer, nrowRes, ncol);

        memcpy(dataPtr, resBuffer, resBufferSize);
        dataPtr = dataPtr + nrowRes;
    }

    delete[] resBuffer;
    delete o;

    DenseMatrix resT(resultData, ncolRes, nrowRes, false);
    return resT.transpose();
}

DenseMatrix* DenseMatrix::copy() {
    float* data = this->_data;
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;

    size_t size = ncol * nrow;

    float* newData = new float[size];
    memcpy(newData, data, size * sizeof(float));

    return new DenseMatrix(newData, nrow, ncol, false);
}


DenseVector* gaussJordanEliminationVector(DenseMatrix *A, DenseVector *vector) {
    // https://martin-thoma.com/solving-linear-equations-with-gaussian-elimination/
    size_t n = A->numRows();
    assert(n == vector->size());
    assert(A->numRows() == A->numCols());

    DenseMatrix *U = A->copy();
    DenseVector *b = vector->copy();

    for (size_t i = 0; i < n - 1; i++) {
        // 1. search for the max value in this col
        float maxel = abs(U->get(i, i));
        size_t maxrow = i;

        for (size_t k = i + 1; k < n; k++) {
            float el = abs(U->get(k, i));
            if (el > maxel) {
                maxel = el;
                maxrow = k;
            }
        }

        // 2. swap the rows
        if (maxrow != i) {
            U->swapRows(i, maxrow);
            b->swap(i, maxrow);
        }

        maxel = U->get(i, i);
        if (abs(maxel) <= 1e-6) {
            throw invalid_argument("the matrix is singular");
        }

        // 3. adjust values according to maxel
        for (size_t k = i + 1; k < n; k++) {
            float m = U->get(k, i) / maxel;
            float currentValue, rowValue;

            U->set(k, i, 0);
            for (size_t j = i + 1; j < n; j++) {
                currentValue = U->get(k, j);
                rowValue = U->get(i, j);
                U->set(k, j, currentValue - m * rowValue);
            }

            currentValue = b->get(k);
            rowValue = b->get(i);
            b->set(k, currentValue - m * rowValue);
        }
    }


    // 4. Solve for upper triangular matrix U
    float* x = new float[n];

    for (int i = n - 1; i >= 0; i--) {
        x[i] = b->get(i) / U->get(i, i);
        for (int k = i - 1; k >= 0; k--) {
            float s = U->get(k, i) * x[i];
            b->set(k, b->get(k) - s);
        }
    }

    delete U;
    delete b;

    return new DenseVector(x, n, false);
}


DenseMatrix* gaussJordanEliminationMatrix(DenseMatrix *A, DenseMatrix *matrix) {
    size_t n = A->numRows();
    assert(A->numRows() == A->numCols());
    assert(n == matrix->numRows());

    size_t bNumCols = matrix->numCols();

    DenseMatrix *U = A->copy();
    DenseMatrix *B = matrix->copy();

    for (size_t i = 0; i < n - 1; i++) {
        // 1. search for the max value in this col
        float maxel = abs(U->get(i, i));
        size_t maxrow = i;

        for (size_t k = i + 1; k < n; k++) {
            float el = abs(U->get(k, i));
            if (el > maxel) {
                maxel = el;
                maxrow = k;
            }
        }

        // 2. swap the rows
        if (maxrow != i) {
            U->swapRows(i, maxrow);
            B->swapRows(i, maxrow);
        }

        maxel = U->get(i, i);
        if (abs(maxel) <= 1e-6) {
            throw invalid_argument("the matrix is singular");
        }

        // 3. adjust values according to maxel
        maxel = U->get(i, i);
        for (size_t k = i + 1; k < n; k++) {
            float rowMultiplier = U->get(k, i) / maxel;
            float currentValue, rowValue;

            U->set(k, i, 0);
            for (size_t j = i + 1; j < n; j++) {
                currentValue = U->get(k, j);
                rowValue = U->get(i, j);
                U->set(k, j, currentValue - rowMultiplier * rowValue);
            }

            for (size_t j = 0; j < bNumCols; j++) {
                currentValue = B->get(k, j);
                rowValue = B->get(i, j);
                B->set(k, j, currentValue - rowMultiplier * rowValue);
            }
            
        }
    }

    // 4. Solve for upper triangular matrix U
    DenseMatrix *X = new DenseMatrix(n, bNumCols);

    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < bNumCols; j++) {
            float newVal = B->get(i, j) / U->get(i, i);
            X->set(i, j, newVal);

            for (int k = i - 1; k >= 0; k--) {
                float s = U->get(k, i) * newVal;
                B->set(k, j, B->get(k, j) - s);
            }
        }
    }

    delete U;
    delete B;

    return X;
}

DenseVector* DenseMatrix::solve(DenseVector *vector) {
    return gaussJordanEliminationVector(this, vector);
}

DenseMatrix* DenseMatrix::solveMatrix(DenseMatrix *matrix) {
    return gaussJordanEliminationMatrix(this, matrix);
}


// LUDecomposition::~LUDecomposition() {
//     cout << "destructor called" << endl;
// }

LUDecomposition DenseMatrix::lu() {
    // https://www.quantstart.com/articles/LU-Decomposition-in-Python-and-NumPy
    size_t n = this->_nrow;
    assert(this->_nrow == this->_ncol);

    DenseMatrix *PA = this->copy();
    DenseMatrix *P = new DenseMatrix(n, n);
    DenseMatrix *L = new DenseMatrix(n, n);
    DenseMatrix *U = new DenseMatrix(n, n);

    // build P: rearrange elements such that max is located on the diagonal
    for (size_t i = 0; i < n; i++) {
        P->set(i, i, 1);
    }

    for (size_t i = 0; i < n; i++) {
        float maxel = abs(PA->get(i, i));
        size_t maxrow = i;

        for (size_t k = i + 1; k < n; k++) {
            float el = abs(PA->get(k, i));
            if (el > maxel) {
                maxel = el;
                maxrow = k;
            }
        }

        if (maxrow != i) {
            PA->swapRows(i, maxrow);
            P->swapRows(i, maxrow);
        }
    }


    for (size_t i = 0; i < n; i++) {
        L->set(i, i, 1);

        for (size_t k = 0; k <= i; k++) {
            float s1 = 0.0;

            for (size_t j = 0; j < k; j++) {
                s1 = s1 + U->get(j, i) * L->get(k, j);
            }

            float p = PA->get(k, i);
            U->set(k, i, p - s1);
        }

        for (size_t k = i + 1; k < n; k++) {
            float s2 = 0.0;

            for (size_t j = 0; j < i; j++) {
                s2 = s2 + U->get(j, i) * L->get(k, j);
            }

            float p = PA->get(k, i);
            float u = U->get(i, i);
            L->set(k, i, (p - s2) / u);
        }
    }

    delete PA;

    LUDecomposition result;
    result.P = P;
    result.L = L;
    result.U = U;
    return result;
}

DenseMatrix* DenseMatrix::inverseGaussJordan() {
    size_t n = this->_nrow;
    assert(this->_nrow == this->_ncol);
    DenseMatrix *eye = DenseMatrix::eye(n);
    DenseMatrix *inv = this->solveMatrix(eye);
    delete eye;
    return inv;
}


DenseMatrix* upperTriangularSolveMatrix(DenseMatrix *U, DenseMatrix *matrix, 
        bool solveInplace) {
    size_t n = U->numRows();
    assert(U->numRows() == U->numCols());
    assert(n == matrix->numRows());
    // check that U is upper triangular ?

    size_t bNumCols = matrix->numCols();
    DenseMatrix *B = solveInplace ? matrix : matrix->copy();
    DenseMatrix *X = new DenseMatrix(n, bNumCols);

    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < bNumCols; j++) {
            float newVal = B->get(i, j) / U->get(i, i);
            X->set(i, j, newVal);

            for (int k = i - 1; k >= 0; k--) {
                float s = U->get(k, i) * newVal;
                B->set(k, j, B->get(k, j) - s);
            }
        }
    }

    if (!solveInplace) {
        delete B;    
    }

    return X;
}

DenseMatrix* lowerTriangularSolveMatrix(DenseMatrix *L, DenseMatrix *matrix, 
        bool solveInplace) {
    size_t n = L->numRows();
    assert(L->numRows() == L->numCols());
    assert(n == matrix->numRows());
    // check that L is lower triangular ?

    size_t bNumCols = matrix->numCols();
    DenseMatrix *B = solveInplace ? matrix : matrix->copy();

    DenseMatrix *X = new DenseMatrix(n, bNumCols);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < bNumCols; j++) {
            float newVal = B->get(i, j) / L->get(i, i);
            X->set(i, j, newVal);

            for (size_t k = i + 1; k < n; k++) {
                float s = L->get(k, i) * newVal;
                B->set(k, j, B->get(k, j) - s);
            }
        }
    }

    if (!solveInplace) {
        delete B;    
    }

    return X;
}

DenseMatrix* DenseMatrix::inverseLU() {
    assert(this->_nrow == this->_ncol);
    LUDecomposition PLU = this->lu();

    DenseMatrix *P = PLU.P, *L = PLU.L, *U = PLU.U;

    DenseMatrix *Y = lowerTriangularSolveMatrix(L, P, true);
    DenseMatrix *X = upperTriangularSolveMatrix(U, Y, true);

    delete P;
    delete L;
    delete U;
    delete Y;

    return X;
}

DenseMatrix* DenseMatrix::inverse() {
    return this->inverseQR();
}

DenseMatrix* DenseMatrix::orthonormalize() {
    assert(this->_nrow == this->_ncol);
    size_t n = this->_nrow;

    DenseMatrix *QT = this->transpose();

    float* data = QT->_data;
    size_t ncol = QT->_ncol;

    arrayUnitize(data, ncol);

    for (size_t i = 1; i < n; i++) {
        float* q_i = &data[i * ncol];

        for (size_t j = 0; j < i; j++) {
            float* q_j = &data[j * ncol];
            float t = -arrayDot(q_i, q_j, ncol);
            addScaledVector(q_i, q_j, t, ncol);
        }

        arrayUnitize(q_i, ncol);        
    }

    DenseMatrix *Q = QT->transpose();
    delete QT;
    return Q;
}


QRDecomposition DenseMatrix::qr() {
    assert(this->_nrow == this->_ncol);
    size_t n = this->_nrow;

    DenseMatrix *QT = this->transpose();
    DenseMatrix *R = new DenseMatrix(n, n);

    float* data = QT->_data;
    size_t ncol = QT->_ncol;

    for (int i = 0; i < n; i++) {
        float* q_i = &data[i * ncol];

        for (int j = 0; j < i; j++) {
            float* q_j = &data[j * ncol];
            float t = arrayDot(q_i, q_j, ncol);
            R->set(j, i, t);
            addScaledVector(q_i, q_j, -t, ncol);
        }

        float norm2 = arrayNorm2(q_i, ncol);
        float norm = sqrt(norm2);
        R->set(i, i, norm);
        arrayScale(q_i, 1 / norm, ncol);    
    }

    QRDecomposition result;
    result.QT = QT;
    result.R = R;
    return result;
}

DenseMatrix* DenseMatrix::inverseQR() {
    assert(this->_nrow == this->_ncol);
    size_t n = this->_nrow;

    QRDecomposition QR = this->qr();

    DenseMatrix *QT = QR.QT, *R = QR.R;
    DenseMatrix *X = upperTriangularSolveMatrix(R, QT, true);

    delete QT;
    delete R;

    return X;
}

size_t DenseMatrix::numRows() {
    return this->_nrow;    
}

size_t DenseMatrix::numCols() {
    return this->_ncol;
}

void DenseMatrix::printMatrix() {    
    size_t nrow = this->_nrow;
    size_t ncol = this->_ncol;

    cout << "nrow: " << nrow << ", ncol: " << ncol << endl;
    for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            cout << this->get(i, j) << " ";
        }
        cout << endl;
    }
}
