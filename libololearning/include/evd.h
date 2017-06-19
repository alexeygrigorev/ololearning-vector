#pragma once

#include "DenseVector.h"
#include "DenseMatrix.h"

namespace olomath {

struct EVD {
    DenseMatrix *V;
    DenseVector *S;
};


struct SVD {
    DenseMatrix *U;
    DenseVector *S;
    DenseMatrix *V;
};

EVD truncatedEVD(DenseMatrix A, int k, float kappa);

SVD truncatedSVD(DenseMatrix A, int k, float kappa);


}