#include "DenseVector.h"
#include <iostream>
#include <cmath>
#include <random>

#include "DenseMatrix.h"
#include "gtest/gtest.h"

using namespace std;

class VectorTestSuite : public ::testing::Test {};


TEST_F(VectorTestSuite, VectorInitializedWithZeros) {
    DenseVector v(15);
    EXPECT_EQ(v.get(1), 0.0);
}

TEST_F(VectorTestSuite, VectorSetAndGet) {
    DenseVector v(15);
    v.set(1, 10.0);
    EXPECT_EQ(v.get(1), 10.0);
}

TEST_F(VectorTestSuite, VectorInitializedWithArray) {
    const size_t size = 4;
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);

    EXPECT_EQ(v.get(0), 0.0);
    EXPECT_EQ(v.get(1), 1.0);
    EXPECT_EQ(v.get(2), 2.0);
    EXPECT_EQ(v.get(3), 3.0);
}

TEST_F(VectorTestSuite, VectorSwap) {
    const size_t size = 4;
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    v.swap(1, 2);

    EXPECT_EQ(v.get(0), 0.0);
    EXPECT_EQ(v.get(1), 2.0);
    EXPECT_EQ(v.get(2), 1.0);
    EXPECT_EQ(v.get(3), 3.0);
}

TEST_F(VectorTestSuite, VectorNorm2) {
    const size_t size = 4;
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    float norm2 = v.norm2();
    float expected = 1 + 4 + 9;

    EXPECT_EQ(norm2, expected);
}

TEST_F(VectorTestSuite, VectorScale) {
    const size_t size = 4;
    float* data = new float[size] { 0, 1, 2, 3 };

    DenseVector v(data, size);
    DenseVector *u = v.scale(2, false);

    float* dataScaled = new float[size] { 0, 2, 4, 6 };
    DenseVector scaled(dataScaled, size);

    float dist = u->distance2(&scaled);
    ASSERT_TRUE(dist <= 1e-6);
}


TEST_F(VectorTestSuite, VectorUnitNormalize) {
    const size_t size = 4;
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    DenseVector *u = v.unitize(false);

    EXPECT_EQ(1 + 4 + 9, v.norm2());   
    ASSERT_TRUE(abs(u->norm2() - 1) < 1e-6);
}


TEST_F(VectorTestSuite, VectorUnitNormalizeInplace) {
    const size_t size = 4;
    float* data = new float[size] { 0.0, 1.0, 2.0, 3.0 };

    DenseVector v(data, size);
    v.unitize(true);

    float norm2 = v.norm2();
    ASSERT_TRUE(abs(norm2 - 1) < 1e-6);
}


TEST_F(VectorTestSuite, VectorDistance2) {
    const size_t size = 4;
    float* data1 = new float[size] { 0.0, 1.0, 2.0, 3.0 };
    float* data2 = new float[size] { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(data1, size), v2(data2, size);
    float dist2 = v1.distance2(&v2);
    float expected = 9 + 4 + 1 + 0;

    EXPECT_EQ(dist2, expected);
}


TEST_F(VectorTestSuite, VectorDot) {
    const size_t size = 4;
    float* data1 = new float[size] { 0.0, 1.0, 2.0, 3.0 };
    float* data2 = new float[size] { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(data1, size), v2(data2, size);
    float dot = v1.dot(&v2);
    float expected = 3 + 6 + 9;

    EXPECT_EQ(dot, expected);
}


TEST_F(VectorTestSuite, VectorSubtract) {
    const size_t size = 4;
    float* data1 = new float[size] { 0.0, 1.0, 2.0, 3.0 };
    float* data2 = new float[size] { 3.0, 3.0, 3.0, 3.0 };

    DenseVector v1(data1, size), v2(data2, size);
    DenseVector *s = v1.subtract(&v2, false);

    EXPECT_EQ(s->get(0), -3);
    EXPECT_EQ(s->get(1), -2);
    EXPECT_EQ(s->get(2), -1);
    EXPECT_EQ(s->get(3), 0);
}