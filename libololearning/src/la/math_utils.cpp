#include <cstddef>
#include <cstring>
#include <cmath>

static float* zeroVector(size_t size) {
    float* data = new float[size];
    memset(data, 0, size * sizeof(float));
    return data;
}

static float* zeroMatrix(size_t nrow, size_t ncol) {
    return zeroVector(nrow * ncol);
}


static float arrayNorm2(float* data, size_t size) {
    float norm2 = 0.0;
    for (size_t i = 0; i < size; i++) {
        float el = data[i];
        norm2 = norm2 + el * el;
    }

    return norm2;
}

static float* copyOrSame(float* array, size_t size, bool inplace) {
    if (inplace == true) {
        return array;
    }

    float* newData = new float[size];
    memcpy(newData, array, size * sizeof(float));
    return newData;
}

static void arrayScale(float *data, float scale, size_t size) {
    for (size_t i = 0; i < size; i++) {
         data[i] = data[i] * scale;
    }
}

static void arrayUnitize(float *data, size_t size) {
    float norm2 = arrayNorm2(data, size);
    float norm = sqrt(norm2);
    arrayScale(data, 1 / norm, size);
}

static float arrayDot(float *u, float *v, size_t size) {
    float dot = 0;
    for (size_t i = 0; i < size; i++) {
        dot = dot + u[i] * v[i];
    }
    return dot;
}

static void addScaledVector(float* target, float *source, float scale, size_t size) {
    for (size_t i = 0; i < size; i++) {
        target[i] = target[i] + scale * source[i];
    }
}
