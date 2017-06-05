/**
@file Tridiagonalization.cpp
@author t-sakai
@date 2017/06/02 create
*/
#include "lmath.h"
#include "Vector.h"
#include "Matrix.h"
#include "Householder.h"

namespace lmath
{
namespace
{
    void updateWithHouseholder(Matrix& m, Vector& v, s32 size)
    {
        s32 offset = m.cols() - size;
        Vector g(size);

        for(s32 i=offset; i<m.cols(); ++i){
            g[i-offset] = 0.0;
            for(s32 j=offset; j<m.cols(); ++j){
                g[i-offset] += m(j,i)*v[j-offset];
            }
        }
        element_type d = dot(g, v, size);
        for(s32 i=0; i<size; ++i){
            g[i] = 2.0f*(g[i] - v[i]*d);
        }
        for(s32 i=offset; i<m.cols(); ++i){
            for(s32 j=offset; j<m.cols(); ++j){
                m(j,i) -= (v[i-offset]*g[j-offset] + g[i-offset]*v[j-offset]);
            }
        }
    }
}

    void tridiagonalize(Matrix& m)
    {
        Vector v(m.cols());
        s32 cols = m.cols()-2;
        s32 size = v.size();
        for(s32 i=0; i<cols; ++i){
            --size;
            for(s32 j=i+1; j<m.cols(); ++j){
                v[j-i-1] = m(j,i);
            }
            element_type norm = householder(v, size);
            if(-LMATH_EPSILON<norm && norm<LMATH_EPSILON){
                continue;
            }
            updateWithHouseholder(m, v, size);
            m(i+1, i) = m(i,i+1) = norm;
            for(s32 j=i+2; j<m.cols(); ++j){
                m(j,i) = m(i,j) = 0.0;
            }
        }
    }

    void tridiagonalize(s32 n, MatrixView& m, element_type* d, element_type* e)
    {
        s32 n2 = n-2;
        s32 n1 = n-1;
        for(s32 k=0; k<n2; ++k){
            element_type* v = &m(0,k);
            d[k] = v[k];
            e[k] = householder(&v[k+1], n1-k);
            if(-LMATH_EPSILON<e[k] && e[k]<LMATH_EPSILON){
                continue;
            }
            for(s32 i=k+1; i<n; ++i){
                element_type sum = 0;
                for(s32 j=k+1; j<i; ++j){
                    sum += m(i,j)*v[j];
                }
                for(s32 j=i; j<n; ++j){
                    sum += m(j,i)*v[j];
                }
                d[i] = sum;
            }
            element_type t = innerproduct(n1-k, &v[k+1], &d[k+1])*0.5;
            for(s32 i=n1; k<i; --i){
                element_type p = v[i];
                element_type q = d[i]-t*p;
                d[i] = q;
                for(s32 j=i; j<n; ++j){
                    m(j,i) -= p*d[j] + q*v[j];
                }
            }
        }//for(s32 k=0

        if(2<=n){
            d[n2] = m(n2,n2);
            e[n2] = m(n1,n2);
        }
        if(1<=n){
            d[n1] = m(n1, n1);
        }
        for(s32 k=n1; 0<=k; --k){
            element_type* v = &m(0,k);
            if(k<n2){
                for(s32 i=k+1; i<n; ++i){
                    element_type* w = &m(0,i);
                    element_type t = innerproduct(n1-k, &v[k+1], &w[k+1]);
                    for(s32 j=k+1; j<n; ++j){
                        w[j] -= t*v[j];
                    }
                }
            }
            for(s32 i=0; i<n; ++i){
                v[i] = 0.0;
            }
            v[k] = 1.0;
        }
    }
}
