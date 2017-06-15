#include "DenseMatrix.h"
#include "testsuite.h"

#include <iostream>
#include <cmath>
#include <random>

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
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);

    EXPECT_EQ(v.get(0), 0.0);
    EXPECT_EQ(v.get(1), 1.0);
    EXPECT_EQ(v.get(2), 2.0);
    EXPECT_EQ(v.get(3), 3.0);
}

TEST_F(TestSuite, VectorSwap) {
    const size_t size = 4;
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    v.swap(1, 2);

    EXPECT_EQ(v.get(0), 0.0);
    EXPECT_EQ(v.get(1), 2.0);
    EXPECT_EQ(v.get(2), 1.0);
    EXPECT_EQ(v.get(3), 3.0);
}

TEST_F(TestSuite, VectorNorm2) {
    const size_t size = 4;
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    float norm2 = v.norm2();
    float expected = 1 + 4 + 9;

    EXPECT_EQ(norm2, expected);
}

TEST_F(TestSuite, VectorScale) {
    const size_t size = 4;
    float* data = new float[size] { 0, 1, 2, 3 };

    DenseVector v(data, size);
    DenseVector u = v.scale(2, false);

    float* dataScaled = new float[size] { 0, 2, 4, 6 };
    DenseVector scaled(dataScaled, size);

    float dist = u.distance2(scaled);
    ASSERT_TRUE(dist <= 1e-6);
}


TEST_F(TestSuite, VectorUnitNormalize) {
    const size_t size = 4;
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    DenseVector u = v.unitize(false);

    EXPECT_EQ(1 + 4 + 9, v.norm2());   
    ASSERT_TRUE(abs(u.norm2() - 1) < 1e-6);
}


TEST_F(TestSuite, VectorUnitNormalizeInplace) {
    const size_t size = 4;
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    v.unitize(true);

    float norm2 = v.norm2();
    ASSERT_TRUE(abs(norm2 - 1) < 1e-6);
}




TEST_F(TestSuite, VectorDistance2) {
    const size_t size = 4;
    float* data1 = new float[size] { 0.0, 1.0, 2.0, 3.0 };
    float* data2 = new float[size] { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(data1, size), v2(data2, size);
    float dist2 = v1.distance2(v2);
    float expected = 9 + 4 + 1 + 0;

    EXPECT_EQ(dist2, expected);
}


TEST_F(TestSuite, VectorDot) {
    const size_t size = 4;
    float* data1 = new float[size] { 0.0, 1.0, 2.0, 3.0 };
    float* data2 = new float[size] { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(data1, size), v2(data2, size);
    float dot = v1.dot(v2);
    float expected = 3 + 6 + 9;

    EXPECT_EQ(dot, expected);
}


TEST_F(TestSuite, VectorSubtract) {
    const size_t size = 4;
    float* data1 = new float[size] { 0.0, 1.0, 2.0, 3.0 };
    float* data2 = new float[size] { 3.0, 3.0, 3.0, 3.0 };

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
    float data[nrow][ncol] = {
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

TEST_F(TestSuite, MatrixClone) {
    float dataA[2][2] = {
        { 1, 2 },
        { 2, 2 },
    };
    DenseMatrix A(&dataA[0][0], 2, 2);
    DenseMatrix copy = A.copy();
    copy.set(0, 0, 2);

    EXPECT_EQ(copy.get(0, 0), 2);
    EXPECT_EQ(A.get(0, 0), 1);
}


TEST_F(TestSuite, MatrixGetColumn) {
    const int nrow = 4, ncol = 3;
    float data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m(&data[0][0], nrow, ncol);
    DenseVector col = m.getColumn(1);

    float colData[nrow] = { 1.0, 4.0, 7.0, 0.0 };
    DenseVector expected(&colData[0], nrow);

    float diff = col.distance2(expected);
    EXPECT_EQ(0, diff);
}

TEST_F(TestSuite, MatrixGetRow) {
    const int nrow = 4, ncol = 3;
    float data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m(&data[0][0], nrow, ncol);
    DenseVector row = m.getRow(1);

    float rowData[ncol] = { 3.0, 4.0, 5.0 };
    DenseVector expected(&rowData[0], ncol);

    float diff = row.distance2(expected);
    EXPECT_EQ(0, diff);
}

TEST_F(TestSuite, MatrixSwapRow) {
    const int nrow = 4, ncol = 3;
    float data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix m(&data[0][0], nrow, ncol);
    m.swapRows(1, 2);

    DenseVector row1 = m.getRow(1);
    float rowData1[ncol] = { 6.0, 7.0, 8.0 };
    DenseVector expected1(&rowData1[0], ncol);
    EXPECT_EQ(0, row1.distance2(expected1));

    DenseVector row2 = m.getRow(2);
    float rowData2[ncol] = { 3.0, 4.0, 5.0 };
    DenseVector expected2(&rowData2[0], ncol);
    EXPECT_EQ(0, row2.distance2(expected2));
}

TEST_F(TestSuite, MatrixSubtract) {
    const int nrow = 2, ncol = 3;
    float data1[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
    };

    float data2[nrow][ncol] = {
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
    float data[nrow][ncol] = {
        {0.0, 1.0},
        {2.0, 3.0},
    };

    DenseMatrix m(&data[0][0], nrow, ncol);
    float norm = m.norm2();
    float expected = 1.0 + 4.0 + 9.0;
    EXPECT_EQ(norm, expected);
}

TEST_F(TestSuite, MatrixDistance2) {
    const int nrow = 2, ncol = 2;
    float data1[nrow][ncol] = {
        {0.0, 1.0},
        {2.0, 3.0},
    };
    float data2[nrow][ncol] = {
        {2.0, 2.0},
        {2.0, 2.0},
    };

    DenseMatrix m1(&data1[0][0], nrow, ncol);
    DenseMatrix m2(&data2[0][0], nrow, ncol);

    float diff = m1.distance2(m2);
    float expected = 4.0 + 1.0 + 1.0;

    EXPECT_EQ(diff, expected);
}


TEST_F(TestSuite, MatrixTranspose) {
    const int nrow = 4, ncol = 3;
    float data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };
    float expectedData[ncol][nrow] = {
        { 0.0, 3.0, 6.0, 9.0 },
        { 1.0, 4.0, 7.0, 0.0 },
        { 2.0, 5.0, 8.0, 1.0 }
    };

    DenseMatrix m(&data[0][0], nrow, ncol);
    DenseMatrix t = m.transpose();

    DenseMatrix expT(&expectedData[0][0], ncol, nrow);

    float dist = t.distance2(expT);
    EXPECT_EQ(dist, 0);
}

TEST_F(TestSuite, MatrixVectorMult) {
    const int nrow = 4, ncol = 3;
    float data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };
    float vdata[ncol] = { 1.0, 2.0, 3.0 };
    float edata[nrow] = { 8., 26., 44., 12. };

    DenseMatrix m(&data[0][0], nrow, ncol);
    DenseVector v(&vdata[0], ncol);

    DenseVector res = m.vmult(v);
    DenseVector expected(&edata[0], nrow);

    float dist = res.distance2(expected);
    EXPECT_EQ(dist, 0);
}


TEST_F(TestSuite, MatrixMatrixMult) {
    float dataA[3][2] = {
        {0.0, 1.0},
        {3.0, 4.0},
        {6.0, 7.0},
    };
    float dataB[2][4] = {
        {0.0, 1.0, 2.0, 3.0},
        {3.0, 4.0, 1.0, 2.0},
    };
    float dataC[3][4] = { 
        {  3.,   4.,   1.,   2.},
        { 12.,  19.,  10.,  17.},
        { 21.,  34.,  19.,  32.}
    };

    DenseMatrix A(&dataA[0][0], 3, 2);
    DenseMatrix B(&dataB[0][0], 2, 4);
    DenseMatrix expC(&dataC[0][0], 3, 4);

    DenseMatrix C = A.mmult(B);

    float dist = C.distance2(expC);
    EXPECT_EQ(dist, 0);
}


TEST_F(TestSuite, MatrixSolve3x3Simple) {
    float dataA[3][3] = {
        { 1, 1, 1 },
        { 0, 1, 1 },
        { 0, 0, 1 }
    };
    float dataB[3] = { 6, 3, 1 };
    float dataXExp[3] = { 3,  2, 1 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseVector b(&dataB[0], 3);
    DenseVector x = A.solve(b);

    DenseVector expX(&dataXExp[0], 3);
    
    float dist = x.distance2(expX);
    ASSERT_TRUE(dist <= 0.005);
}



TEST_F(TestSuite, MatrixSolve3x3) {
    float dataA[3][3] = {
        { 1, 2, 1 },
        { 3, 8, 1 },
        { 9, 4, 1 }
    };
    float dataB[3] = { 2, 12, 2 };
    float dataXExp[3] = { -0.454,  1.818, -1.181 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseVector b(&dataB[0], 3);
    DenseVector x = A.solve(b);
    
    DenseVector expX(&dataXExp[0], 3);

    float dist = x.distance2(expX);
    ASSERT_TRUE(dist <= 0.005);
}

TEST_F(TestSuite, MatrixSolve3x3Matrix) {
    float dataA[3][3] = {
        { 1, 2, 1 },
        { 3, 8, 1 },
        { 9, 4, 1 }
    };
    float dataB[3][3] = {
        { 2,  2,  2  },
        { 12, 12, 12 },
        { 2,  2,  2  }
    };
    float dataXExp[3][3] = {
        { -0.454, -0.454, -0.454 },
        {  1.818,  1.818,  1.818 }, 
        { -1.181, -1.181, -1.181 }
    };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseMatrix B(&dataB[0][0], 3, 3);
    DenseMatrix X = A.solveMatrix(B);

    DenseMatrix expX(&dataXExp[0][0], 3, 3);

    float dist = X.distance2(expX);
    ASSERT_TRUE(dist <= 0.005);
}

TEST_F(TestSuite, MatrixSolve3x3Matrix_SmallB) {
    float dataA[3][3] = {
        { 1, 2, 1 },
        { 3, 8, 1 },
        { 9, 4, 1 }
    };
    float dataB[3][1] = {
        { 2  },
        { 12 },
        { 2  }
    };
    float dataXExp[3][1] = {
        { -0.454 },
        {  1.818 },
        { -1.181 }
    };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseMatrix B(&dataB[0][0], 3, 1);
    DenseMatrix X = A.solveMatrix(B);
    cout << X.numRows() << ", " << X.numCols() << endl;

    DenseMatrix expX(&dataXExp[0][0], 3, 1);

    float dist = X.distance2(expX);
    ASSERT_TRUE(dist <= 0.005);
}

TEST_F(TestSuite, MatrixSolve3x3Matrix_WideB) {
    float dataA[3][3] = {
        { 1, 2, 1 },
        { 3, 8, 1 },
        { 9, 4, 1 }
    };
    float dataB[3][6] = {
        { 2,  2,  2,  1, 0, 0 },
        { 12, 12, 12, 0, 1, 0 },
        { 2,  2,  2,  0, 0, 1 }
    };
    float dataXExp[3][6] = {
        { -0.454, -0.454, -0.454, -0.090, -0.045,  0.136 },
        {  1.818,  1.818,  1.818, -0.136,  0.181, -0.045 }, 
        { -1.181, -1.181, -1.181,  1.363, -0.313, -0.045 }
    };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseMatrix B(&dataB[0][0], 3, 6);
    DenseMatrix X = A.solveMatrix(B);

    DenseMatrix expX(&dataXExp[0][0], 3, 6);

    float dist = X.distance2(expX);
    ASSERT_TRUE(dist <= 0.005);
}


TEST_F(TestSuite, MatrixSolve3x3MatrixInverse) {
    float dataA[3][3] = {
        { 1, 2, 1 },
        { 3, 8, 1 },
        { 9, 4, 1 }
    };
    float dataB[3][3] = {
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 }
    };
    float dataXExp[3][3] = {
        { -0.090, -0.045,  0.136 },
        { -0.136,  0.181, -0.045 }, 
        {  1.363, -0.313, -0.045 }
    };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseMatrix B(&dataB[0][0], 3, 3);
    DenseMatrix X = A.solveMatrix(B);

    DenseMatrix expX(&dataXExp[0][0], 3, 3);

    float dist = X.distance2(expX);
    ASSERT_TRUE(dist <= 0.005);
}

TEST_F(TestSuite, MatrixInverse2x2) {
    float dataA[2][2] = {
        { 1, 2 },
        { 2, 2 },
    };
    float dataAInv[2][2] = {
        { -1.0,  1.0 },
        {  1.0, -0.5 },
    };

    DenseMatrix A(&dataA[0][0], 2, 2);
    DenseMatrix AInv = A.inverse();

    DenseMatrix expAInv(&dataAInv[0][0], 2, 2);

    float dist = AInv.distance2(expAInv);
    ASSERT_TRUE(dist <= 0.000001);
}

TEST_F(TestSuite, VectorSolve3x3_2) {
    float dataA[3][3] = {
        { 60, 91, 26 },
        { 60,  3, 75 },
        { 45, 90, 31 }
    };
    float dataB[3] = { 1, 0, 0 };
    float dataXExp[3] = { 0.053, -0.012, -0.042 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseVector b(&dataB[0], 3);
    DenseVector x = A.solve(b);
    
    DenseVector expX(&dataXExp[0], 3);

    float dist = x.distance2(expX);
    ASSERT_TRUE(dist <= 0.001);
}

TEST_F(TestSuite, VectorSolve3x3_Singular) {
    float dataA[3][3] = {
        { 91, 91, 26 },
        {  3,  3, 75 },
        { 90, 90, 31 }
    };
    float dataB[3] = { 1, 0, 0 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseVector b(&dataB[0], 3);

    try {
        DenseVector x = A.solve(b);
        FAIL();    
    } catch (const invalid_argument& ex) {
        // pass
    }
}


TEST_F(TestSuite, MatrixInverse3x3) {
    float dataA[3][3] = {
        { 60, 91, 26 },
        { 60,  3, 75 },
        { 45, 90, 31 }
    };
    float dataAInv[3][3] = {
        {  0.053,  0.003, -0.054 },
        { -0.012, -0.005,  0.023 },
        { -0.042,  0.010,  0.042 }
    };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseMatrix AInv = A.inverse();

    DenseMatrix expAInv(&dataAInv[0][0], 3, 3);

    float dist = AInv.distance2(expAInv);
    ASSERT_TRUE(dist <= 0.001);
}


TEST_F(TestSuite, MatrixSolve3x3_Singular) {
    float dataA[3][3] = {
        { 91, 91, 26 },
        {  3,  3, 75 },
        { 90, 90, 31 }
    };
    float dataB[3] = { 1, 0, 0 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseMatrix B(&dataB[0], 3, 1);

    try {
        DenseMatrix X = A.solveMatrix(B);
        FAIL();    
    } catch (const invalid_argument& ex) {
        // pass
    }
}

TEST_F(TestSuite, MatrixInverse50x50) {
    int seed = 10;
    std::mt19937_64 random_gen(seed);
    std::uniform_real_distribution<float> unif_random(0, 1);

    int n = 50;
    int size = n * n;
    float* data = new float[size];
    for (int i = 0; i < size; i++) {
        data[i] = unif_random(random_gen);
    }

    DenseMatrix A(data, n, n);

    DenseMatrix AInv = A.inverse();
    DenseMatrix AAinv = A.mmult(AInv);

    DenseMatrix I = DenseMatrix::eye(n);

    float dist = AAinv.distance2(I);
    cout << "error: " << dist << endl;
    ASSERT_TRUE(dist <= 1e-6);
}


TEST_F(TestSuite, MatrixInverseGJ100x100) {
    int seed = 10;
    std::mt19937_64 random_gen(seed);
    std::uniform_real_distribution<float> unif_random(0, 1);

    int n = 100;
    int size = n * n;
    float* data = new float[size];
    for (int i = 0; i < size; i++) {
        data[i] = unif_random(random_gen);
    }

    DenseMatrix A(data, n, n);

    DenseMatrix AInv = A.inverseGaussJordan();
    DenseMatrix AAinv = A.mmult(AInv);

    DenseMatrix I = DenseMatrix::eye(n);

    float dist = AAinv.distance2(I);
    cout << "error: " << dist << endl;
    ASSERT_TRUE(dist <= 1e-6);
}


TEST_F(TestSuite, LUDecomposition4x4) {
    float dataA[4][4] = {
        { 7,  3, -1,  2 },
        { 3,  8,  1, -4 },
        {-1,  1,  4, -1 },
        { 2, -4, -1,  6 }
    };

    float dataP[4][4] = {
        { 1, 0, 0, 0 },
        { 0, 1, 0, 0 },
        { 0, 0, 1, 0 },
        { 0, 0, 0, 1 }
    };
    float dataL[4][4] = {
        { 1.   ,  0.   ,  0.   ,  0. },
        { 0.428,  1.   ,  0.   ,  0. },
        {-0.142,  0.212,  1.   ,  0. },
        { 0.285, -0.723,  0.089,  1. }
    };
    float dataU[4][4] = {
        { 7.,  3.   , -1.   ,  2.    },
        { 0.,  6.714,  1.428, -4.857 },
        { 0.,  0.   ,  3.553,  0.319 },
        { 0.,  0.   ,  0.   ,  1.886 }
    };

    DenseMatrix expP(&dataP[0][0], 4, 4);
    DenseMatrix expL(&dataL[0][0], 4, 4);
    DenseMatrix expU(&dataU[0][0], 4, 4);

    DenseMatrix A(&dataA[0][0], 4, 4);
    LUDecomposition result = A.lu();

    float distP = result.P->distance2(expP);
    ASSERT_TRUE(distP <= 0.00001);

    float distL = result.L->distance2(expL);
    ASSERT_TRUE(distL <= 0.001);

    float distU = result.U->distance2(expU);
    ASSERT_TRUE(distU <= 0.001);
}

TEST_F(TestSuite, MatrixLU50x50) {
    int seed = 10;
    std::mt19937_64 random_gen(seed);
    std::uniform_real_distribution<float> unif_random(0, 1);

    int n = 50;
    int size = n * n;
    float* data = new float[size];
    for (int i = 0; i < size; i++) {
        data[i] = unif_random(random_gen);
    }

    DenseMatrix A(data, n, n);

    LUDecomposition PLU = A.lu();
    DenseMatrix *P = PLU.P, *L = PLU.L, *U = PLU.U;

    DenseMatrix PA = P->mmult(A);
    DenseMatrix LU = L->mmult(*U);

    float dist = PA.distance2(LU);
    cout << "error: " << dist << endl;
    ASSERT_TRUE(dist <= 1e-6);

    delete P;
    delete L;
    delete U;
}


TEST_F(TestSuite, MatrixInverseLU100x100) {
    int seed = 100;
    std::mt19937_64 random_gen(seed);
    std::uniform_real_distribution<float> unif_random(0, 1);

    int n = 100;
    int size = n * n;
    float* data = new float[size];
    for (int i = 0; i < size; i++) {
        data[i] = unif_random(random_gen);
    }

    DenseMatrix A(data, n, n);

    DenseMatrix AInv = A.inverseLU();
    DenseMatrix AAinv = A.mmult(AInv);

    DenseMatrix I = DenseMatrix::eye(n);

    float dist = AAinv.distance2(I);
    cout << "error: " << dist << endl;
    ASSERT_TRUE(dist <= 1e-4);
}


TEST_F(TestSuite, MatrixGramSchmidt) {
    int seed = 100;
    std::mt19937_64 random_gen(seed);
    std::uniform_real_distribution<float> unif_random(0, 1);

    int n = 100;
    int size = n * n;
    float* data = new float[size];
    for (int i = 0; i < size; i++) {
        data[i] = unif_random(random_gen);
    }

    DenseMatrix A(data, n, n);

    DenseMatrix Q = A.orthonormalize();

    DenseMatrix QT = Q.transpose();
    DenseMatrix QTQ = QT.mmult(Q);
    DenseMatrix QQT = Q.mmult(QT);

    DenseMatrix I = DenseMatrix::eye(n);

    float dist1 = QTQ.distance2(I);
    cout << dist1 << endl;
    ASSERT_TRUE(dist1 <= 1e-4);

    float dist2 = QQT.distance2(I);
    ASSERT_TRUE(dist2 <= 1e-4);

}


