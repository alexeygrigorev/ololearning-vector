#include <cstddef>
#include <cstring>
#include <cmath>

#include "math_utils.h"

namespace olomath {


float* zeroVector(size_t size) {
    float* data = new float[size];
    memset(data, 0, size * sizeof(float));
    return data;
}

float* zeroMatrix(size_t nrow, size_t ncol) {
    return zeroVector(nrow * ncol);
}


float arrayNorm2(float* data, size_t size) {
    float norm2 = 0.0;
    for (size_t i = 0; i < size; i++) {
        float el = data[i];
        norm2 = norm2 + el * el;
    }

    return norm2;
}

float* copyOrSame(float* array, size_t size, bool inplace) {
    if (inplace == true) {
        return array;
    }

    float* newData = new float[size];
    memcpy(newData, array, size * sizeof(float));
    return newData;
}

void arrayScale(float *data, float scale, size_t size) {
    for (size_t i = 0; i < size; i++) {
         data[i] = data[i] * scale;
    }
}

void arrayUnitize(float *data, size_t size) {
    float norm2 = arrayNorm2(data, size);
    float norm = sqrt(norm2);
    arrayScale(data, 1 / norm, size);
}

float arrayDot(float *u, float *v, size_t size) {
    float dot = 0;
    for (size_t i = 0; i < size; i++) {
        dot = dot + u[i] * v[i];
    }
    return dot;
}

void matrixVectorProduct(float *M, float *b, float *res, size_t nrow, size_t ncol) {
    for (size_t i = 0; i < nrow; i++) {
        float *row = &M[i * ncol];
        res[i] = arrayDot(row, b, ncol);
    }
}

void addScaledVector(float* target, float *source, float scale, size_t size) {
    for (size_t i = 0; i < size; i++) {
        target[i] = target[i] + scale * source[i];
    }
}


} //namespace
