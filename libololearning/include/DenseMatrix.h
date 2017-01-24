#include <cstddef>

class DenseMatrix {
public: 
    DenseMatrix(size_t nrow, size_t ncol);
    static DenseMatrix fromRowArray(double* data, size_t nrow, size_t ncol);
    static DenseMatrix eye(size_t n);
    static DenseMatrix ones(size_t nrow, size_t ncol);

    DenseMatrix(double* data, size_t nrow, size_t ncol);
    double get(size_t row, size_t col);
    void set(size_t row, size_t col, double val);

    double norm2();

    double* getColumn(size_t col);
    double* getRow(size_t row);

	DenseMatrix subtract(DenseMatrix other, bool inplace); 

    DenseMatrix transpose(); 

    void printMatrix();

private:
	void init(double* data, size_t nrow, size_t ncol);
    size_t nrow;
    size_t ncol;
    double* data;
};