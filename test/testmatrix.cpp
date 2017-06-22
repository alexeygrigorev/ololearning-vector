#include "DenseVector.h"
#include <iostream>
#include <cmath>
#include <random>

#include "DenseMatrix.h"
#include "evd.h"
#include "gtest/gtest.h"

using namespace std;

class MatrixTestSuite : public ::testing::Test {};


TEST_F(MatrixTestSuite, MatrixInitializedWithZeros) {
    DenseMatrix M(10, 15);
    EXPECT_EQ(M.get(0, 3), 0.0);
}

TEST_F(MatrixTestSuite, MatrixSetAndGet) {
    DenseMatrix M(10, 15);
    M.set(3, 3, 5.0);
    EXPECT_EQ(M.get(3, 3), 5.0);
}


TEST_F(MatrixTestSuite, MatrixInitializedWithArray) {
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

TEST_F(MatrixTestSuite, MatrixClone) {
    float dataA[2][2] = {
        { 1, 2 },
        { 2, 2 },
    };
    DenseMatrix A(&dataA[0][0], 2, 2);
    DenseMatrix *copy = A.copy();
    copy->set(0, 0, 2);

    EXPECT_EQ(copy->get(0, 0), 2);
    EXPECT_EQ(A.get(0, 0), 1);
}


TEST_F(MatrixTestSuite, MatrixGetColumn) {
    const int nrow = 4, ncol = 3;
    float data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix M(&data[0][0], nrow, ncol);
    DenseVector *col = M.getColumn(1);

    float colData[nrow] = { 1.0, 4.0, 7.0, 0.0 };
    DenseVector expected(&colData[0], nrow);

    float diff = col->distance2(&expected);
    EXPECT_EQ(0, diff);
}

TEST_F(MatrixTestSuite, MatrixGetRow) {
    const int nrow = 4, ncol = 3;
    float data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix M(&data[0][0], nrow, ncol);
    DenseVector *row = M.getRow(1);

    float rowData[ncol] = { 3.0, 4.0, 5.0 };
    DenseVector expected(&rowData[0], ncol);

    float diff = row->distance2(&expected);
    EXPECT_EQ(0, diff);
}

TEST_F(MatrixTestSuite, MatrixSwapRow) {
    const int nrow = 4, ncol = 3;
    float data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix M(&data[0][0], nrow, ncol);
    M.swapRows(1, 2);


    DenseVector *row1 = M.getRow(1);
    float rowData1[ncol] = { 6.0, 7.0, 8.0 };
    DenseVector expected1(&rowData1[0], ncol);

    float dist1 = row1->distance2(&expected1);
    EXPECT_EQ(0, dist1);


    DenseVector *row2 = M.getRow(2);
    float rowData2[ncol] = { 3.0, 4.0, 5.0 };
    DenseVector expected2(&rowData2[0], ncol);

    float dist2 = row2->distance2(&expected2);
    EXPECT_EQ(0, dist2);
}

TEST_F(MatrixTestSuite, MatrixSubtract) {
    const int nrow = 2, ncol = 3;
    float data1[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
    };

    float data2[nrow][ncol] = {
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };

    DenseMatrix M1(&data1[0][0], nrow, ncol);
    DenseMatrix M2(&data2[0][0], nrow, ncol);

    DenseMatrix *res = M1.subtract(&M2, false);

    EXPECT_EQ(res->get(0, 0), -6.0);
    EXPECT_EQ(res->get(0, 1), -6.0);
    EXPECT_EQ(res->get(0, 2), -6.0);
    EXPECT_EQ(res->get(1, 0), -6.0);
    EXPECT_EQ(res->get(1, 1),  4.0);
    EXPECT_EQ(res->get(1, 2),  4.0);
}

TEST_F(MatrixTestSuite, MatrixNorm2) {
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

TEST_F(MatrixTestSuite, MatrixDistance2) {
    const int nrow = 2, ncol = 2;
    float data1[nrow][ncol] = {
        {0.0, 1.0},
        {2.0, 3.0},
    };
    float data2[nrow][ncol] = {
        {2.0, 2.0},
        {2.0, 2.0},
    };

    DenseMatrix M1(&data1[0][0], nrow, ncol);
    DenseMatrix M2(&data2[0][0], nrow, ncol);

    float diff = M1.distance2(&M2);
    float expected = 4.0 + 1.0 + 1.0;

    EXPECT_EQ(diff, expected);
}


TEST_F(MatrixTestSuite, MatrixTranspose) {
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
    DenseMatrix *t = m.transpose();

    DenseMatrix expT(&expectedData[0][0], ncol, nrow);

    float dist = t->distance2(&expT);
    EXPECT_EQ(dist, 0);
}

TEST_F(MatrixTestSuite, MatrixVectorMult) {
    const int nrow = 4, ncol = 3;
    float data[nrow][ncol] = {
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {9.0, 0.0, 1.0},
    };
    float vdata[ncol] = { 1.0, 2.0, 3.0 };
    float edata[nrow] = { 8., 26., 44., 12. };

    DenseMatrix M(&data[0][0], nrow, ncol);
    DenseVector v(&vdata[0], ncol);

    DenseVector *res = M.vmult(&v);
    DenseVector expected(&edata[0], nrow);

    float dist = res->distance2(&expected);
    EXPECT_EQ(dist, 0);
}


TEST_F(MatrixTestSuite, MatrixMatrixMult) {
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

    DenseMatrix *C = A.mmult(&B);

    float dist = C->distance2(&expC);
    EXPECT_EQ(dist, 0);
}


TEST_F(MatrixTestSuite, MatrixSolve3x3Simple) {
    float dataA[3][3] = {
        { 1, 1, 1 },
        { 0, 1, 1 },
        { 0, 0, 1 }
    };
    float dataB[3] = { 6, 3, 1 };
    float dataXExp[3] = { 3,  2, 1 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseVector b(&dataB[0], 3);
    DenseVector *x = A.solve(&b);

    DenseVector expX(&dataXExp[0], 3);
    
    float dist = x->distance2(&expX);
    ASSERT_TRUE(dist <= 0.005);
}



TEST_F(MatrixTestSuite, MatrixSolve3x3) {
    float dataA[3][3] = {
        { 1, 2, 1 },
        { 3, 8, 1 },
        { 9, 4, 1 }
    };
    float dataB[3] = { 2, 12, 2 };
    float dataXExp[3] = { -0.454,  1.818, -1.181 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseVector b(&dataB[0], 3);
    DenseVector *x = A.solve(&b);
    
    DenseVector expX(&dataXExp[0], 3);

    float dist = x->distance2(&expX);
    ASSERT_TRUE(dist <= 0.005);
}

TEST_F(MatrixTestSuite, MatrixSolve3x3Matrix) {
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
    DenseMatrix *X = A.solveMatrix(&B);

    DenseMatrix expX(&dataXExp[0][0], 3, 3);

    float dist = X->distance2(&expX);
    ASSERT_TRUE(dist <= 0.005);
}

TEST_F(MatrixTestSuite, MatrixSolve3x3Matrix_SmallB) {
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
    DenseMatrix *X = A.solveMatrix(&B);

    DenseMatrix expX(&dataXExp[0][0], 3, 1);

    float dist = X->distance2(&expX);
    ASSERT_TRUE(dist <= 0.005);
}

TEST_F(MatrixTestSuite, MatrixSolve3x3Matrix_WideB) {
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
    DenseMatrix *X = A.solveMatrix(&B);

    DenseMatrix expX(&dataXExp[0][0], 3, 6);

    float dist = X->distance2(&expX);
    ASSERT_TRUE(dist <= 0.005);
}


TEST_F(MatrixTestSuite, MatrixSolve3x3MatrixInverse) {
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
    DenseMatrix *X = A.solveMatrix(&B);

    DenseMatrix expX(&dataXExp[0][0], 3, 3);

    float dist = X->distance2(&expX);
    ASSERT_TRUE(dist <= 0.005);
}

TEST_F(MatrixTestSuite, MatrixInverse2x2) {
    float dataA[2][2] = {
        { 1, 2 },
        { 2, 2 },
    };
    float dataAInv[2][2] = {
        { -1.0,  1.0 },
        {  1.0, -0.5 },
    };

    DenseMatrix A(&dataA[0][0], 2, 2);
    DenseMatrix *AInv = A.inverse();

    DenseMatrix expAInv(&dataAInv[0][0], 2, 2);

    float dist = AInv->distance2(&expAInv);
    ASSERT_TRUE(dist <= 0.000001);
}

TEST_F(MatrixTestSuite, VectorSolve3x3_2) {
    float dataA[3][3] = {
        { 60, 91, 26 },
        { 60,  3, 75 },
        { 45, 90, 31 }
    };
    float dataB[3] = { 1, 0, 0 };
    float dataXExp[3] = { 0.053, -0.012, -0.042 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseVector b(&dataB[0], 3);
    DenseVector *x = A.solve(&b);
    
    DenseVector expX(&dataXExp[0], 3);

    float dist = x->distance2(&expX);
    ASSERT_TRUE(dist <= 0.001);
}

TEST_F(MatrixTestSuite, VectorSolve3x3_Singular) {
    float dataA[3][3] = {
        { 91, 91, 26 },
        {  3,  3, 75 },
        { 90, 90, 31 }
    };
    float dataB[3] = { 1, 0, 0 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseVector b(&dataB[0], 3);

    try {
        A.solve(&b);
        FAIL();    
    } catch (const invalid_argument& ex) {
        // pass
    }
}


TEST_F(MatrixTestSuite, MatrixInverse3x3) {
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
    DenseMatrix *AInv = A.inverse();

    DenseMatrix expAInv(&dataAInv[0][0], 3, 3);

    float dist = AInv->distance2(&expAInv);
    ASSERT_TRUE(dist <= 0.001);
}


TEST_F(MatrixTestSuite, MatrixSolve3x3_Singular) {
    float dataA[3][3] = {
        { 91, 91, 26 },
        {  3,  3, 75 },
        { 90, 90, 31 }
    };
    float dataB[3] = { 1, 0, 0 };

    DenseMatrix A(&dataA[0][0], 3, 3);
    DenseMatrix B(&dataB[0], 3, 1);

    try {
        DenseMatrix *X = A.solveMatrix(&B);
        FAIL();    
    } catch (const invalid_argument& ex) {
        // pass
    }
}

TEST_F(MatrixTestSuite, MatrixInverse50x50) {
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

    DenseMatrix *AInv = A.inverse();
    DenseMatrix *AAinv = A.mmult(AInv);

    DenseMatrix *I = DenseMatrix::eye(n);

    float dist = AAinv->distance2(I);
    ASSERT_TRUE(dist <= 1e-6);
}


TEST_F(MatrixTestSuite, MatrixInverseGJ100x100) {
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

    DenseMatrix *AInv = A.inverseGaussJordan();
    DenseMatrix *AAinv = A.mmult(AInv);

    DenseMatrix *I = DenseMatrix::eye(n);

    float dist = AAinv->distance2(I);
    ASSERT_TRUE(dist <= 1e-6);
}


TEST_F(MatrixTestSuite, LUDecomposition4x4) {
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

    float distP = result.P->distance2(&expP);
    ASSERT_TRUE(distP <= 0.00001);

    float distL = result.L->distance2(&expL);
    ASSERT_TRUE(distL <= 0.001);

    float distU = result.U->distance2(&expU);
    ASSERT_TRUE(distU <= 0.001);
}

TEST_F(MatrixTestSuite, MatrixLU50x50) {
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

    DenseMatrix *PA = P->mmult(&A);
    DenseMatrix *LU = L->mmult(U);

    float dist = PA->distance2(LU);
    ASSERT_TRUE(dist <= 1e-6);

    delete P;
    delete L;
    delete U;
}


TEST_F(MatrixTestSuite, MatrixInverseLU100x100) {
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

    DenseMatrix *AInv = A.inverseLU();
    DenseMatrix *AAinv = A.mmult(AInv);

    DenseMatrix *I = DenseMatrix::eye(n);

    float dist = AAinv->distance2(I);
    ASSERT_TRUE(dist <= 1e-4);
}


TEST_F(MatrixTestSuite, MatrixGramSchmidt) {
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

    DenseMatrix *Q = A.orthonormalize();

    DenseMatrix *QT = Q->transpose();
    DenseMatrix *QTQ = QT->mmult(Q);
    DenseMatrix *QQT = Q->mmult(QT);

    DenseMatrix *I = DenseMatrix::eye(n);

    float dist1 = QTQ->distance2(I);
    ASSERT_TRUE(dist1 <= 1e-4);

    float dist2 = QQT->distance2(I);
    ASSERT_TRUE(dist2 <= 1e-4);

}


TEST_F(MatrixTestSuite, MatrixQR_100x100) {
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

    QRDecomposition qr = A.qr();
    DenseMatrix *QT = qr.QT, *R = qr.R;

    DenseMatrix *Q = QT->transpose();
    DenseMatrix *QR = Q->mmult(R);
    float dist0 = QR->distance2(&A);
    ASSERT_TRUE(dist0 <= 1e-4);

    DenseMatrix *QTQ = QT->mmult(Q);
    DenseMatrix *QQT = Q->mmult(QT);
    DenseMatrix *I = DenseMatrix::eye(n);

    float dist1 = QTQ->distance2(I);
    ASSERT_TRUE(dist1 <= 1e-4);

    float dist2 = QQT->distance2(I);
    ASSERT_TRUE(dist2 <= 1e-4);

    DenseMatrix *QTA = QT->mmult(&A);
    float dist3 = QTA->distance2(R);
    ASSERT_TRUE(dist3 <= 1e-4);

}


TEST_F(MatrixTestSuite, MatrixInverseQR100x100) {
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

    DenseMatrix *AInv = A.inverseQR();
    DenseMatrix *AAinv = A.mmult(AInv);

    DenseMatrix *I = DenseMatrix::eye(n);

    float dist = AAinv->distance2(I);
    ASSERT_TRUE(dist <= 1e-4);
}
