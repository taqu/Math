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
    class Vector;
    class VectorView;

    element_type householder(Vector& v);
    element_type householder(Vector& v, s32 size);
    element_type householder(element_type* v, s32 size);
}
#endif //INC_LMATH_HOUSEHOLDER_H__
