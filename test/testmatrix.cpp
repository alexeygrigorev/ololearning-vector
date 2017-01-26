#include "DenseMatrix.h"
#include "testsuite.h"

using namespace std;

TestSuite::TestSuite() {}
TestSuite::~TestSuite() {};

void TestSuite::SetUp() {};
void TestSuite::TearDown() {};


TEST_F(TestSuite, VectorInitializedWithZeros) {
    DenseVector v(15);
    EXPECT_EQ(v.get(1), 0.0);
}

TEST_F(TestSuite, VectorSetAndGet) {
    DenseVector v(15);
    v.set(1, 10.0);
    EXPECT_EQ(v.get(1), 10.0);
}

TEST_F(TestSuite, VectorInitializedWithArray) {
    const size_t size = 4;
    double data[size] = { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(&data[0], size);

    EXPECT_EQ(v.get(0), 0.0);
    EXPECT_EQ(v.get(1), 1.0);
    EXPECT_EQ(v.get(2), 2.0);
    EXPECT_EQ(v.get(3), 3.0);
}

TEST_F(TestSuite, VectorNorm2) {
    const size_t size = 4;
    double data[size] = { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(&data[0], size);
    double norm2 = v.norm2();
    double expected = 1 + 4 + 9;

    EXPECT_EQ(norm2, expected);
}

TEST_F(TestSuite, VectorDistance2) {
    const size_t size = 4;
    double data1[size] = { 0.0, 1.0, 2.0, 3.0 };
    double data2[size] = { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(&data1[0], size), v2(&data2[0], size);
    double dist2 = v1.distance2(v2);
    double expected = 9 + 4 + 1;

    EXPECT_EQ(dist2, expected);
}


TEST_F(TestSuite, VectorDot) {
    const size_t size = 4;
    double data1[size] = { 0.0, 1.0, 2.0, 3.0 };
    double data2[size] = { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(&data1[0], size), v2(&data2[0], size);
    double dot = v1.dot(v2);
    double expected = 3 + 6 + 9;

    EXPECT_EQ(dot, expected);
}


TEST_F(TestSuite, VectorSubtract) {
    const size_t size = 4;
    double data1[size] = { 0.0, 1.0, 2.0, 3.0 };
    double data2[size] = { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(&data1[0], size), v2(&data2[0], size);
    DenseVector s = v1.subtract(v2, false);

    EXPECT_EQ(s.get(0), -3);
    EXPECT_EQ(s.get(1), -2);
    EXPECT_EQ(s.get(2), -1);
    EXPECT_EQ(s.get(3), 0);
}



TEST_F(TestSuite, MatrixInitializedWithZeros) {
    DenseMatrix m(10, 15);
    EXPECT_EQ(m.get(0, 3), 0.0);
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

    DenseMatrix m(&data[0][0], nrow, ncol);
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

    DenseMatrix m(&data[0][0], nrow, ncol);
    DenseVector col = m.getColumn(1);

    double colData[nrow] = { 1.0, 4.0, 7.0, 0.0 };
    DenseVector expected(&colData[0], nrow);

    double diff = col.distance2(expected);
    EXPECT_EQ(0, diff);
}

TEST_F(TestSuite, MatrixGetRow) {
    const int nrow = 4, ncol = 3;
    double data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m(&data[0][0], nrow, ncol);
    DenseVector row = m.getRow(1);

    double rowData[ncol] = { 3.0, 4.0, 5.0 };
    DenseVector expected(&rowData[0], ncol);

    double diff = row.distance2(expected);
    EXPECT_EQ(0, diff);
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

    DenseMatrix m1(&data1[0][0], nrow, ncol);
    DenseMatrix m2(&data2[0][0], nrow, ncol);

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

    DenseMatrix m(&data[0][0], nrow, ncol);
    double norm = m.norm2();
    double expected = 1.0 + 4.0 + 9.0;
    EXPECT_EQ(norm, expected);
}

TEST_F(TestSuite, MatrixDistance2) {
    const int nrow = 2, ncol = 2;
    double data1[nrow][ncol] = {
        {0.0, 1.0},
        {2.0, 3.0},
    };
    double data2[nrow][ncol] = {
        {2.0, 2.0},
        {2.0, 2.0},
    };

    DenseMatrix m1(&data1[0][0], nrow, ncol);
    DenseMatrix m2(&data2[0][0], nrow, ncol);

    double diff = m1.distance2(m2);
    double expected = 4.0 + 1.0 + 1.0;

    EXPECT_EQ(diff, expected);
}


TEST_F(TestSuite, MatrixTranspose) {
    const int nrow = 4, ncol = 3;
    double data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };
    double expectedData[ncol][nrow] = {
        { 0.0, 3.0, 6.0, 9.0 },
        { 1.0, 4.0, 7.0, 0.0 },
        { 2.0, 5.0, 8.0, 1.0 }
    };

    DenseMatrix m(&data[0][0], nrow, ncol);
    DenseMatrix t = m.transpose();

    DenseMatrix expT(&expectedData[0][0], ncol, nrow);

    double dist = t.distance2(expT);
    EXPECT_EQ(dist, 0);
}
