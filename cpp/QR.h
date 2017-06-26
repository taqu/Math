#ifndef INC_LMATH_QR_H__
#define INC_LMATH_QR_H__
/**
@file QR.h
@author t-sakai
@date 2017/06/02 create
*/
#include "lmath.h"

namespace lmath
{
    class VectorView;
    class Matrix;
    class MatrixView;
#if 0
    void qr_algorithm(Matrix& m, element_type epsilon=LMATH_EPSILON);

    void qr_decomp(MatrixView& m);
    void xtoq(MatrixView& q, const MatrixView& x);

    /**
    @return 0Ç»ÇÁê¨å˜
    */
    s32 eigen(s32 n, MatrixView& m, VectorView& d, VectorView& e, f32 epsilon=LMATH_EPSILON, s32 maxIteration=10);
#endif
}
#endif //INC_LMATH_QR_H__
