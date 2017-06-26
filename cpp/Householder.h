#ifndef INC_LMATH_HOUSEHOLDER_H__
#define INC_LMATH_HOUSEHOLDER_H__
/**
@file Householder.h
@author t-sakai
@date 2017/06/02 create
*/
#include "lmath.h"

namespace lmath
{
    class VectorView;
    class VectorStepView;
    class MatrixView;

    element_type householder(element_type* v, s32 size);

    /**
    @warning v‚Ì“à—e‚ð”j‰ó‚·‚é
    */
    void householder_matrix(MatrixView& m, element_type* v, s32 k, s32 size);

    /**
    @warning v‚Ì“à—e‚ð”j‰ó‚·‚é
    */
    void householder_matrix(MatrixView& m, VectorView& v, s32 k);

    /**
    @warning v‚Ì“à—e‚ð”j‰ó‚·‚é
    */
    void householder_matrix(MatrixView& m, VectorStepView& v, s32 k);
}
#endif //INC_LMATH_HOUSEHOLDER_H__
