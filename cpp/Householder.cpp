/**
@file Householder.cpp
@author t-sakai
@date 2017/06/02 create
*/
#include <math.h>
#include "Householder.h"
#include "Matrix.h"

namespace lmath
{
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

    void householder_matrix(MatrixView& m, element_type* v, s32 k, s32 size)
    {
        LASSERT(0<size);
        LASSERT(0<=k);
        LASSERT(size==m.rows());
        LASSERT(size==m.cols());

        s32 d = size-k;
        v += k;
        element_type norm = lmath::sqrt(innerproduct(d,v,v));

        if(norm<LMATH_EPSILON){
            m.identity();
            return;
        }
        element_type alpha = 0.0<=v[0]? -norm : norm;

        v[0] = lmath::sqrt(0.5*(1.0-v[0]/alpha));
        element_type weight = 1.0/(-2.0*alpha*v[0]);
        for(s32 i=1; i<d; ++i){
            v[i] *= weight;
        }

        for(s32 i=0; i<m.rows(); ++i){
            for(s32 j=0; j<k; ++j){
                m(i,j)=0.0;
                m(j,j) = 1.0;
                for(s32 l=k; l<m.cols(); ++l){
                    m(j,l) = 0.0;
                }
            }
        }

        for(s32 i=k; i<m.rows(); ++i){
            for(s32 j=k; j<m.cols(); ++j){
                m(i,j) = -2.0*v[i-k]*v[j-k];
            }
            m(i,i) += 1.0;
        }
    }

    void householder_matrix(MatrixView& m, VectorView& v, s32 k)
    {
        householder_matrix(m, &v[0], k, v.size());
    }

    void householder_matrix(MatrixView& m, VectorStepView& v, s32 k)
    {
        LASSERT(0<v.size());
        LASSERT(0<=k && k<v.size());
        LASSERT(v.size()==m.rows());
        LASSERT(v.size()==m.cols());

        s32 d = v.size()-k;
        element_type norm = lmath::sqrt(innerproduct(k,v,v));
        element_type alpha = 0.0<=v[k]? -norm : norm;

        if(norm<LMATH_EPSILON){
            m.identity();
            return;
        }
        v[k] = lmath::sqrt(0.5*(1.0-v[k]/alpha));
        element_type weight = 1.0/(-2.0*alpha*v[k]);
        for(s32 i=1; i<d; ++i){
            v[k+i] *= weight;
        }

        for(s32 i=0; i<m.rows(); ++i){
            for(s32 j=0; j<k; ++j){
                m(i,j)=0.0;
            }
        }
        for(s32 i=0; i<k; ++i){
            m(i,i) = 1.0;
            for(s32 j=k; j<m.cols(); ++j){
                m(i,j) = 0.0;
            }
        }

        for(s32 i=k; i<m.rows(); ++i){
            for(s32 j=k; j<m.cols(); ++j){
                m(i,j) = -2.0*v[i-k]*v[j-k];
            }
            m(i,i) += 1.0;
        }
    }
}
