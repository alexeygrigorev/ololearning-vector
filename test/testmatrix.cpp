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
    double* data = new double[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);

    EXPECT_EQ(v.get(0), 0.0);
    EXPECT_EQ(v.get(1), 1.0);
    EXPECT_EQ(v.get(2), 2.0);
    EXPECT_EQ(v.get(3), 3.0);
}

TEST_F(TestSuite, VectorSwap) {
    const size_t size = 4;
    double* data = new double[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    v.swap(1, 2);

    EXPECT_EQ(v.get(0), 0.0);
    EXPECT_EQ(v.get(1), 2.0);
    EXPECT_EQ(v.get(2), 1.0);
    EXPECT_EQ(v.get(3), 3.0);
}

TEST_F(TestSuite, VectorNorm2) {
    const size_t size = 4;
    double* data = new double[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    double norm2 = v.norm2();
    double expected = 1 + 4 + 9;

    EXPECT_EQ(norm2, expected);
}

TEST_F(TestSuite, VectorDistance2) {
    const size_t size = 4;
    double* data1 = new double[size] { 0.0, 1.0, 2.0, 3.0 };
    double* data2 = new double[size] { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(data1, size), v2(data2, size);
    double dist2 = v1.distance2(v2);
    double expected = 9 + 4 + 1 + 0;

    EXPECT_EQ(dist2, expected);
}


TEST_F(TestSuite, VectorDot) {
    const size_t size = 4;
    double* data1 = new double[size] { 0.0, 1.0, 2.0, 3.0 };
    double* data2 = new double[size] { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(data1, size), v2(data2, size);
    double dot = v1.dot(v2);
    double expected = 3 + 6 + 9;

    EXPECT_EQ(dot, expected);
}


TEST_F(TestSuite, VectorSubtract) {
    const size_t size = 4;
    double* data1 = new double[size] { 0.0, 1.0, 2.0, 3.0 };
    double* data2 = new double[size] { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(data1, size), v2(data2, size);
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

TEST_F(TestSuite, MatrixSwapRow) {
    const int nrow = 4, ncol = 3;
    double data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m(&data[0][0], nrow, ncol);
    m.swapRows(1, 2);

    DenseVector row1 = m.getRow(1);
    double rowData1[ncol] = { 6.0, 7.0, 8.0 };
    DenseVector expected1(&rowData1[0], ncol);
    EXPECT_EQ(0, row1.distance2(expected1));

    DenseVector row2 = m.getRow(2);
    double rowData2[ncol] = { 3.0, 4.0, 5.0 };
    DenseVector expected2(&rowData2[0], ncol);
    EXPECT_EQ(0, row2.distance2(expected2));
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

TEST_F(TestSuite, MatrixVectorMult) {
    const int nrow = 4, ncol = 3;
    double data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };
    double vdata[ncol] = { 1.0, 2.0, 3.0 };
    double edata[nrow] = { 8., 26., 44., 12. };

    DenseMatrix m(&data[0][0], nrow, ncol);
    DenseVector v(&vdata[0], ncol);

    DenseVector res = m.vmult(v);
    DenseVector expected(&edata[0], nrow);

    double dist = res.distance2(expected);
    EXPECT_EQ(dist, 0);
}


TEST_F(TestSuite, MatrixMatrixMult) {
    double dataA[3][2] = {
        {0.0, 1.0},
        {3.0, 4.0},
        {6.0, 7.0},
    };
    double dataB[2][4] = {
        {0.0, 1.0, 2.0, 3.0},
        {3.0, 4.0, 1.0, 2.0},
    };
    double dataC[3][4] = { 
        {  3.,   4.,   1.,   2.},
        { 12.,  19.,  10.,  17.},
        { 21.,  34.,  19.,  32.}
    };

    DenseMatrix A(&dataA[0][0], 3, 2);
    DenseMatrix B(&dataB[0][0], 2, 4);
    DenseMatrix expC(&dataC[0][0], 3, 4);

    DenseMatrix C = A.mmult(B);

    double dist = C.distance2(expC);
    EXPECT_EQ(dist, 0);
}

TEST_F(TestSuite, MatrixSolve3x3) {
    double dataA[3][3] = {
        { 1, 2, 1 },
        { 3, 8, 1 },
        { 9, 4, 1 }
    };
    double dataB[3] = { 2, 12, 2 };
    double dataXExp[3] = { -0.454,  1.818, -1.181 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseVector b(&dataB[0], 3);
    DenseVector x = A.solve(b);
    
    DenseVector expX(&dataXExp[0], 3);
    
    cout << "x:" << endl;
    x.printVector();

    cout << "expected x:" << endl;
    expX.printVector();
    double dist = x.distance2(expX);
    cout << dist << endl;

    ASSERT_TRUE(dist <= 0.005);
}



// TEST_F(TestSuite, MatrixInverse2x2) {
//     double dataA[2][2] = {
//         { 1, 2 },
//         { 2, 2 },
//     };
//     double dataAInv[2][2] = {
//         { -1.0,  1.0 },
//         {  1.0, -0.5 },
//     };

//     DenseMatrix A(&dataA[0][0], 2, 2);
//     DenseMatrix AInv = A.inverse();
//     AInv.printMatrix();

//     DenseMatrix expAInv(&dataAInv[0][0], 2, 2);

//     double dist = AInv.distance2(expAInv);
//     cout << dist << endl;
//     ASSERT_TRUE(dist <= 0.000001);

//}

TEST_F(TestSuite, MatrixClone) {
    double dataA[2][2] = {
        { 1, 2 },
        { 2, 2 },
    };
    DenseMatrix A(&dataA[0][0], 2, 2);
    DenseMatrix copy = A.copy();
    copy.set(0, 0, 2);

    EXPECT_EQ(copy.get(0, 0), 2);
    EXPECT_EQ(A.get(0, 0), 1);
}

