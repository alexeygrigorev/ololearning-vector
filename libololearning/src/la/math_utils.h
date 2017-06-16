
namespace olomath {

float* zeroVector(size_t size);
float* zeroMatrix(size_t nrow, size_t ncol);
float arrayNorm2(float* data, size_t size);
float* copyOrSame(float* array, size_t size, bool inplace);
void arrayScale(float *data, float scale, size_t size);
void arrayUnitize(float *data, size_t size);
float arrayDot(float *u, float *v, size_t size);
void addScaledVector(float* target, float *source, float scale, size_t size);

}