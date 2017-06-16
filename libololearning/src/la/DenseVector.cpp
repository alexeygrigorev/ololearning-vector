#include <cstddef>
#include <cstring>
#include <iostream>
#include <cmath>

#include <stdexcept>
#include <cassert>

#include "DenseVector.h"
#include "math_utils.h"

using namespace std;
using namespace olomath;


DenseVector::DenseVector(size_t size) {
    this->_data = zeroVector(size);
    this->_size = size;
}

DenseVector::DenseVector(float* data, size_t size) {
    this->_data = data;
    this->_size = size;
}

DenseVector::~DenseVector() {
    // delete[] this->_data;
}

float DenseVector::get(size_t i) {
    size_t size = this->_size;
    assert(i < size);
    float* data = this->_data;
    return data[i];
}

void DenseVector::set(size_t i, float val) {
    size_t size = this->_size;
    assert(i < size);
    float* data = this->_data;
    data[i] = val;
}

void DenseVector::swap(size_t i, size_t j) {
    size_t size = this->_size;
    assert(i < size);
    assert(j < size);
    float* data = this->_data;

    float tmp = data[i];
    data[i] = data[j];
    data[j] = tmp;
}


float DenseVector::norm2() {
    size_t size = this->_size;
    float* data = this->_data;
    return arrayNorm2(data, size);
}

float DenseVector::distance2(DenseVector other) {
    size_t size = this->_size;
    assert(size == other._size);

    float* data = this->_data;
    float* oData = other._data;

    float total = 0.0;
    for (size_t i = 0; i < size; i++) {
        float d = data[i] - oData[i];
        total = total + d * d;
    }

    return total;
}

float DenseVector::dot(DenseVector other) {
    size_t size = this->_size;
    assert(size == other._size);

    float* data = this->_data;
    float* oData = other._data;

    return arrayDot(data, oData, size);
}

DenseVector DenseVector::scale(float scalar, bool inplace) {
    size_t size = this->_size;
    float* data = copyOrSame(this->_data, size, inplace);
    arrayScale(data, scalar, size);

    if (inplace == true) {
        return *this;
    } else {
        return DenseVector(data, size);
    }    
}

DenseVector DenseVector::project(DenseVector u) {
    float vu = this->dot(u);
    float uNorm2 = u.norm2();
    float x = vu / uNorm2;
    return u.scale(x, false);
}

DenseVector DenseVector::unitize(bool inplace) {
    float normInv = 1 / sqrt(this->norm2());
    return this->scale(normInv, inplace);
}

DenseVector DenseVector::subtract(DenseVector other, bool inplace) {
    size_t size = this->_size;
    assert(size == other._size);

    float* newData = copyOrSame(this->_data, size, inplace);
    float* data = this->_data;
    float* oData = other._data;

    for (size_t i = 0; i < size; i++) {
        newData[i] = data[i] - oData[i];
    }

    if (inplace == true) {
        return *this;
    } else {
        return DenseVector(newData, size);
    }
}



size_t DenseVector::size() {
    return this->_size;
}

float* DenseVector::getData() {
    return this->_data;
}

DenseVector DenseVector::copy() {
    float* data = this->_data;
    size_t size = this->_size;

    float* newData = new float[size];
    memcpy(newData, data, size * sizeof(float));

    return DenseVector(newData, size);
}

void DenseVector::printVector() {    
    size_t size = this->_size;
    float* data = this->_data;

    cout << "size: " << size << endl;
    for (size_t i = 0; i < size; i++) {
        cout << data[i] << " ";
    }

    cout << endl;
}
