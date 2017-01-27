#include <cstddef>
#include <cstring>
#include <iostream>

#include <cassert>

#include "DenseMatrix.h"

using namespace std;

double* zeroVector(size_t size) {
    double* data = new double[size];
    memset(data, 0, size * sizeof(data));
    return data;
}

double* zeroMatrix(size_t nrow, size_t ncol) {
    return zeroVector(nrow * ncol);
}

DenseVector::DenseVector(size_t size) {
    this->_data = zeroVector(size);
    this->_size = size;
}

DenseVector::DenseVector(double* data, size_t size) {
    this->_data = data;
    this->_size = size;
}

DenseVector::~DenseVector() {
    // delete[] this->_data;
}

double DenseVector::get(size_t i) {
    size_t size = this->_size;
    assert(i < size);
    double* data = this->_data;
    return data[i];
}

void DenseVector::set(size_t i, double val) {
    size_t size = this->_size;
    assert(i < size);
    double* data = this->_data;
    data[i] = val;
}

double DenseVector::norm2() {
    size_t size = this->_size;

    double* data = this->_data;
    double norm2 = 0.0;
    for (size_t i = 0; i < size; i++) {
        double el = data[i];
        norm2 = norm2 + el * el;
    }

    return norm2;
}


double DenseVector::distance2(DenseVector other) {
    size_t size = this->_size;
    assert(size == other._size);

    double* data = this->_data;
    double* oData = other._data;

    double total = 0.0;
    for (size_t i = 0; i < size; i++) {
        double d = data[i] - oData[i];
        total = total + d * d;
    }

    return total;
}


double DenseVector::dot(DenseVector other) {
    size_t size = this->_size;
    assert(size == other._size);

    double* data = this->_data;
    double* oData = other._data;

    double total = 0.0;
    for (size_t i = 0; i < size; i++) {
        total = total + data[i] * oData[i];
    }

    return total;
}

double* copyOrSame(double* array, size_t size, bool inplace) {
    if (inplace == true) {
        return array;
    }

    double* newData = new double[size];
    memcpy(newData, array, size * sizeof(double*));
    return newData;
}

DenseVector DenseVector::subtract(DenseVector other, bool inplace) {
    size_t size = this->_size;
    assert(size == other._size);

    double* newData = copyOrSame(this->_data, size, inplace);
    double* data = this->_data;
    double* oData = other._data;

    for (size_t i = 0; i < size; i++) {
        newData[i] = data[i] - oData[i];
    }

    if (inplace == true) {
        return *this;
    } else {
        DenseVector result(newData, size);
        return result;
    }
}

size_t DenseVector::size() {
    return this->_size;
}

double* DenseVector::getData() {
    return this->_data;
}

void DenseMatrix::init(double* data, size_t nrow, size_t ncol) {
    this->_nrow = nrow;
    this->_ncol = ncol;
    this->_data = data;
}

DenseMatrix::DenseMatrix(double* data, size_t nrow, size_t ncol) {
    this->init(data, nrow, ncol);
}

DenseMatrix::DenseMatrix(size_t nrow, size_t ncol) {
    double* data = zeroMatrix(nrow, ncol);
    this->init(data, nrow, ncol);
}

DenseMatrix::~DenseMatrix() {
    // delete[] this->_data;
}


DenseMatrix DenseMatrix::eye(size_t n) {
    DenseMatrix m(n, n);

    for (size_t i = 0; i < n; i++) {
        m.set(i, i, 0);
    }

    return m;
}


double DenseMatrix::get(size_t row, size_t col) {
    size_t ncol = this->_ncol;
    return this->_data[row * ncol + col];
}

void DenseMatrix::set(size_t row, size_t col, double val) {
    size_t ncol = this->_ncol;
    double* data = this->_data;
    data[row * ncol + col] = val;
}

DenseVector DenseMatrix::getColumn(size_t col) {
    size_t nrow = this->_nrow;

    double* colData = new double[nrow];

    for (size_t r = 0; r < nrow; r++) {
        colData[r] = this->get(r, col);
    }

    return DenseVector(colData, nrow);
}

DenseVector DenseMatrix::getRow(size_t row) {
    double* data = this->_data;
    size_t ncol = this->_ncol;
    double* rowData = &data[row * ncol];
    return DenseVector(rowData, ncol);
}


DenseMatrix DenseMatrix::subtract(DenseMatrix other, bool inplace) {
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;
    assert(ncol == other._ncol);
    assert(nrow == other._nrow);

    size_t size = ncol * nrow;
    double* newData = copyOrSame(this->_data, size, inplace);
    double* data = this->_data;
    double* oData = other._data;

    for (size_t i = 0; i < size; i++) {
        newData[i] = data[i] - oData[i];
    }

    if (inplace == true) {
        return *this;
    } else {
        DenseMatrix result(newData, nrow, ncol);
        return result;
    }
}

double DenseMatrix::norm2() {
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;
    size_t size = ncol * nrow;

    double* data = this->_data;
    double norm2 = 0.0;
    for (size_t i = 0; i < size; i++) {
        double el = data[i];
        norm2 = norm2 + el * el;
    }

    return norm2;
}

double DenseMatrix::distance2(DenseMatrix other) {
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;
    assert(ncol == other._ncol);
    assert(nrow == other._nrow);

    size_t size = ncol * nrow;

    double* data = this->_data;
    double* oData = other._data;

    double total = 0.0;
    for (size_t i = 0; i < size; i++) {
        double d = data[i] - oData[i];
        total = total + d * d;
    }

    return total;
}

DenseMatrix DenseMatrix::transpose() {
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;

    DenseMatrix t(ncol, nrow);

    for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            double val = this->get(i, j);
            t.set(j, i, val);
        }
    }

    return t;
}


DenseVector DenseMatrix::vmult(DenseVector vec) {
    size_t nrow = this->_nrow;
    assert(this->_ncol == vec.size());

    double* result = new double[nrow];

    for (size_t i = 0; i < nrow; i++) {
        DenseVector row = this->getRow(i);
        result[i] = row.dot(vec);
    }

    return DenseVector(result, nrow);
}

DenseMatrix DenseMatrix::mmult(DenseMatrix other) {
    assert(this->_ncol == other._nrow);
    size_t nrowRes = this->_nrow;
    size_t ncolRes = other._ncol;

    size_t size = nrowRes * ncolRes;
    double* resultData = new double[size];
    double* dataPtr = resultData;

    DenseMatrix o = other.transpose();

    for (size_t i = 0; i < ncolRes; i++) {
        DenseVector oCol = o.getRow(i);
        DenseVector colRes = this->vmult(oCol);
        memcpy(dataPtr, colRes.getData(), nrowRes * sizeof(double));

        dataPtr = dataPtr + nrowRes;
    }

    DenseMatrix resT(resultData, ncolRes, nrowRes);
    DenseMatrix result = resT.transpose();
    return result;
}

DenseVector DenseMatrix::solve(DenseVector b) {
    // Gaussian elimination
    return b;
}

DenseMatrix::LUDecomposition DenseMatrix::lu() {
    DenseMatrix::LUDecomposition result;
    result.L = this;
    result.U = this;
    return result;
}

DenseMatrix DenseMatrix::inverse() {
    return *this;
}

DenseMatrix DenseMatrix::clone() {
    double* data = this->_data;
    size_t ncol = this->_ncol;
    size_t nrow = this->_nrow;

    size_t size = ncol * nrow;

    double* newData = new double[size];
    memcpy(newData, data, size * sizeof(double));

    return DenseMatrix(newData, nrow, ncol);
}

void DenseVector::printVector() {    
    size_t size = this->_size;
    double* data = this->_data;

    cout << "size: " << size << endl;
    for (size_t i = 0; i < size; i++) {
        cout << data[i] << " ";
    }

    cout << endl;
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
