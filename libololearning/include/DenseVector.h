#include <cstddef>

class DenseVector {
public:
    DenseVector(size_t size);
    DenseVector(float* data, size_t size);
    DenseVector(float* data, size_t size, bool external);

    ~DenseVector();

    float get(size_t i);
    void set(size_t i, float val);
    void swap(size_t i, size_t j);

    float norm2();
    DenseVector* unitize(bool inplace);

    float distance2(DenseVector *other);
    float dot(DenseVector *other);

    DenseVector* scale(float scalar, bool inplace);
    DenseVector* subtract(DenseVector *other, bool inplace);

    DenseVector* project(DenseVector *u);

    void printVector();

    size_t size();
    float* getData();

    DenseVector* copy();
private:
    size_t _size;
    float* _data;
    bool _external;
};
