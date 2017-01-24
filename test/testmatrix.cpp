#include "DenseMatrix.h"
#include "testsuite.h"

using namespace std;

TestSuite::TestSuite() {}
TestSuite::~TestSuite() {};

void TestSuite::SetUp() {};
void TestSuite::TearDown() {};

TEST_F(TestSuite, MatrixInitializedWithZeros) {
    DenseMatrix m(10, 15);
    EXPECT_EQ(m.get(0, 0), 0.0);
}

TEST_F(TestSuite, MatrixSetAndGet) {
    DenseMatrix m(10, 15);
    m.set(3, 3, 5.0);
    EXPECT_EQ(m.get(3, 3), 5.0);
}


TEST_F(TestSuite, MatrixInitializedWithArray) {
    const int nrow = 4, ncol = 3;
    double data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m = DenseMatrix::fromRowArray(&data[0][0], nrow, ncol);
    m.printMatrix();

    EXPECT_EQ(m.get(0, 0), 0.0);
    EXPECT_EQ(m.get(1, 2), 5.0);
    EXPECT_EQ(m.get(3, 0), 9.0);
    EXPECT_EQ(m.get(3, 2), 1.0);
}

TEST_F(TestSuite, MatrixGetColumn) {
    const int nrow = 4, ncol = 3;
    double data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m = DenseMatrix::fromRowArray(&data[0][0], nrow, ncol);
    double* col = m.getColumn(1);

    EXPECT_EQ(col[0], 1.0);
    EXPECT_EQ(col[1], 4.0);
    EXPECT_EQ(col[2], 7.0);
    EXPECT_EQ(col[3], 0.0);
}

TEST_F(TestSuite, MatrixGetRow) {
    const int nrow = 4, ncol = 3;
    double data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m = DenseMatrix::fromRowArray(&data[0][0], nrow, ncol);
    double* row = m.getRow(1);

    EXPECT_EQ(row[0], 3.0);
    EXPECT_EQ(row[1], 4.0);
    EXPECT_EQ(row[2], 5.0);
}

TEST_F(TestSuite, MatrixSubtract) {
    const int nrow = 2, ncol = 3;
    double data1[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
    };

    double data2[nrow][ncol] = {
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m1 = DenseMatrix::fromRowArray(&data1[0][0], nrow, ncol);
    DenseMatrix m2 = DenseMatrix::fromRowArray(&data2[0][0], nrow, ncol);

    DenseMatrix res = m1.subtract(m2, false);

    EXPECT_EQ(res.get(0, 0), -6.0);
    EXPECT_EQ(res.get(0, 1), -6.0);
    EXPECT_EQ(res.get(0, 2), -6.0);
    EXPECT_EQ(res.get(1, 0), -6.0);
    EXPECT_EQ(res.get(1, 1), 4.0);
    EXPECT_EQ(res.get(1, 2), 4.0);
}

TEST_F(TestSuite, MatrixNorm2) {
    const int nrow = 2, ncol = 2;
    double data[nrow][ncol] = {
        {0.0, 1.0},
        {2.0, 3.0},
    };

    DenseMatrix m = DenseMatrix::fromRowArray(&data[0][0], nrow, ncol);
    double norm = m.norm2();
    double expected = 1.0 + 4.0 + 9.0;
    EXPECT_EQ(norm, expected);
}

/*
TEST_F(TestSuite, MatrixTranspose) {
    const int nrow = 4, ncol = 3;
    double data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m = DenseMatrix::fromRowArray(&data[0][0], nrow, ncol);

    EXPECT_EQ(row[0], 3.0);
    EXPECT_EQ(row[1], 4.0);
    EXPECT_EQ(row[2], 5.0);
}
*/