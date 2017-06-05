#ifndef INC_LMATH_TRIDIAGONALIZE_H__
#define INC_LMATH_TRIDIAGONALIZE_H__
/**
@file Tridiagonalize.h
@author t-sakai
@date 2017/06/02 create
*/
#include "lmath.h"

namespace lmath
{
    class VectorView;
    class Matrix;
    class MatrixView;

    void tridiagonalize(Matrix& m);
    void tridiagonalize(s32 n, MatrixView& m, element_type* d, element_type* e);
}
#endif //INC_LMATH_TRIDIAGONALIZE_H__
