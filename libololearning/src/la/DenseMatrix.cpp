#include <cstddef>
#include <cstring>
#include <iostream>

#include <cassert>

#include "DenseMatrix.h"

using namespace std;


void DenseMatrix::init(double* data, size_t nrow, size_t ncol) {
    this->nrow = nrow;
    this->ncol = ncol;
    this->data = data;
}

DenseMatrix::DenseMatrix(double* data, size_t nrow, size_t ncol) {
    this->init(data, nrow, ncol);
}

double* zeros(size_t nrow, size_t ncol) {
    double* data = new double[nrow * ncol];
    memset(data, 0, sizeof(data));
    return data;
}

DenseMatrix::DenseMatrix(size_t nrow, size_t ncol) {
    double* data = zeros(nrow, ncol);
    this->init(data, nrow, ncol);
}

DenseMatrix DenseMatrix::fromRowArray(double* data, size_t nrow, size_t ncol) {
    double* result = new double[nrow * ncol];

    for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            result[i + nrow * j] = data[i * ncol + j];
        }
    }

    DenseMatrix m(result, nrow, ncol);
    return m;
}

DenseMatrix DenseMatrix::eye(size_t n) {
    DenseMatrix m(n, n);

    for (size_t i = 0; i < n; i++) {
        m.set(i, i, 0);
    }

    return m;
}


double DenseMatrix::get(size_t row, size_t col) {
    size_t ncol = this->ncol;
    return this->data[row + col * nrow];
}

void DenseMatrix::set(size_t row, size_t col, double val) {
    size_t ncol = this->ncol;
    double* data = this->data;
    data[row + col * nrow] = val;
}

double* DenseMatrix::getColumn(size_t col) {
    double* data = this->data;
    size_t nrow = this->nrow;
    return &data[col * nrow];
}

double* DenseMatrix::getRow(size_t row) {
    size_t ncol = this->ncol;

    double* rowData = new double[ncol];

    for (size_t c = 0; c < ncol; c++) {
        rowData[c] = this->get(row, c);
    }

    return rowData;
}


DenseMatrix DenseMatrix::subtract(DenseMatrix other, bool inplace) {
    size_t ncol = this->ncol;
    size_t nrow = this->nrow;
    assert(ncol == other.ncol);
    assert(nrow == other.nrow);

    size_t size = ncol * nrow;
    double* data = this->data;
    double* oData = other.data;

    double* newData;
    if (inplace == true) {
        newData = data;
    } else {
        newData = new double[size];
        memcpy(newData, data, size * sizeof(double));
    }

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
    size_t ncol = this->ncol;
    size_t nrow = this->nrow;
    size_t size = ncol * nrow;

    double* data = this->data;
    double norm2 = 0.0;
    for (size_t i = 0; i < size; i++) {
        double el = data[i];
        norm2 = norm2 + el * el;
    }

    return norm2;
}

DenseMatrix DenseMatrix::transpose() {
    size_t ncol = this->ncol;
    size_t nrow = this->nrow;

    DenseMatrix t(ncol, nrow);

    for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            double val = this->get(i, j);
            t.set(j, i, val);
        }
    }

    return t;
}


void DenseMatrix::printMatrix() {    
    size_t nrow = this->nrow;
    size_t ncol = this->ncol;

    cout << "nrow: " << nrow << ", ncol: " << ncol << endl;
    for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            cout << this->get(i, j) << " ";
        }
        cout << endl;
    }
}

