#ifndef INC_LMATH_SVD_H__
#define INC_LMATH_SVD_H__
/**
@file SVD.h
@author t-sakai
@date 2017/06/02 create
*/
#include "lmath.h"
#include "Vector.h"
#include "Matrix.h"

namespace lmath
{
    void bidiagonalization(MatrixView& U, MatrixView& B, MatrixView& V, MatrixView& A);

    void rotGivens(element_type& c, element_type& s, element_type& r, element_type f, element_type g);

    void msweep(MatrixView& B);
    void msweep2(MatrixView& B);
    void vsweep(VectorView& d, VectorView& e);
    /**
    @return ê¨å˜Ç»ÇÁtrue
    */
    bool svd(VectorView& d, VectorView& e, s32& iterations, element_type epsilon=LMATH_EPSILON, s32 maxIterationFactor=500);
 }
#endif //INC_LMATH_SVD_H__
