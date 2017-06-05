/**
@file Householder.cpp
@author t-sakai
@date 2017/06/02 create
*/
#include <math.h>
#include "Householder.h"
#include "Vector.h"

namespace lmath
{
    element_type householder(Vector& v)
    {
        return householder(v, v.size());
    }

    element_type householder(Vector& v, s32 size)
    {
        element_type norm = lmath::sqrt(dot(v,v,size));
        if(LMATH_EPSILON<=norm){
            if(v[0]<0.0f){
                norm = -norm;
            }
            v[0] += norm;
            element_type weight = 1.0/lmath::sqrt(2.0*norm*v[0]);
            for(s32 i=0; i<size; ++i){
                v[i] *= weight;
            }
        }
        return -norm;
    }

    element_type householder(element_type* v, s32 size)
    {
        element_type norm = lmath::sqrt(innerproduct(size, v, v));
        if(LMATH_EPSILON<=norm){
            if(v[0]<0.0f){
                norm = -norm;
            }
            v[0] += norm;
            element_type weight = 1.0/lmath::sqrt(norm*v[0]);
            for(s32 i=0; i<size; ++i){
                v[i] *= weight;
            }
        }
        return -norm;
    }
}
