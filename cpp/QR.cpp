/**
@file QR.cpp
@author t-sakai
@date 2017/06/02 create
*/
#include "QR.h"
#include <math.h>
#include <stdio.h>
#include "Vector.h"
#include "Matrix.h"

namespace lmath
{
#if 0
namespace
{
    element_type calcEigenValue22(element_type a00, element_type a01, element_type a10, element_type a11)
    {
        element_type b = a11 + a00;
        element_type determinant = b*b - 4.0*(a11*a00 - a01*a10);

        element_type eigen0;
        element_type eigen1;
        if(determinant<=0.0){
            eigen0 = b * 0.5;
            eigen1 = b * 0.5;
        }else{
            eigen0 = (b + lmath::sqrt(determinant)) * 0.5;
            eigen1 = (b - lmath::sqrt(determinant)) * 0.5;
        }

        return (absolute(a11 - eigen0) < absolute(a11 - eigen1))
            ? eigen0
            : eigen1;
    }

    void nextAQ(Matrix& a, MatrixView& q, s32 currentSize)
    {
        for(s32 i=0; i<currentSize-1; i++){
            element_type aii = a(i,i);
            element_type aii2 = a(i, i+1);
            element_type alpha = lmath::sqrt(aii*aii + aii2*aii2);
            element_type s,c;
            if(0.0<alpha){
                element_type ialpha = 1.0/alpha;
                s = a(i,i+1)*ialpha;
                c = a(i,i)*ialpha;
            }else{
                s = c = alpha;
            }

            for(s32 j=i+1; j<currentSize; ++j){
                element_type tmp = -a(j,i)*s + a(j,i+1)*c;
                a(j,i) = a(j,i)*c + a(j, i+1)*s;
                a(j,i+1) = tmp;
            }

            for(s32 j=0; j<currentSize; ++j){
                element_type tmp = -q(i,j)*s + q(i+1, j)*c;
                q(i,j) = q(i,j)*c + q(i+1,j)*s;
                q(i+1,j) = tmp;
            }

            a(i,i) = alpha;
            a(i,i+1) = 0.0;
        }
    }
}


    void qr_algorithm(Matrix& m, element_type epsilon)
    {
        tridiagonalize(m);

        s32 currentSize = m.cols();
        Matrix q(currentSize, currentSize);
        Vector v(currentSize);
        MatrixView qview(q);

        while(1<currentSize){
            s32 c2 = currentSize-2;
            s32 c1 = currentSize-1;
            if(absolute(m(c2, c1))<epsilon){
                --currentSize;
                continue;
            }

            element_type u = calcEigenValue22(m(c2,c2), m(c1,c2), m(c2,c1), m(c1,c1));
            for(s32 i=0; i<currentSize; ++i){
                m(i,i) -= u;
            }
            qview.setSize(currentSize, currentSize);
            qview.setIdentity();
            nextAQ(m, qview, currentSize);

            for(s32 i=0; i<currentSize; ++i){
                for(s32 j=0; j<currentSize; ++j){
                    element_type sum = 0.0;
                    for(s32 k=i; k<currentSize; ++k){
                        sum += m(k,i)*qview(j,k);
                    }
                    v[j] = sum;
                }
                for(s32 j=0; j<currentSize; ++j){
                    m(j,i) = v[j];
                }
            }
            for(s32 i=0; i<currentSize; ++i){
                m(i,i) += u;
            }
        }
    }

    void qr_decomp(MatrixView& m)
    {
        for(s32 i=0; i<m.rows(); ++i){
            VectorView v = m.getRow(i);
            element_type u = lmath::sqrt(innerproduct(m.cols()-i, &v[i], &v[i]));
            if(v[i]<0.0){
                u = -u;
            }
            v[i] += u;
            element_type t = 1.0/(v[i]*u);
            for(s32 j=i+1; j<m.rows(); ++j){
                VectorView w = m.getRow(j);
                element_type s = t * innerproduct(m.cols()-i, &v[i], &w[i]);
                for(s32 k=i; k<m.cols(); ++k){
                    w[k] -= s*v[k];
                }
            }
            v[i] = -u;
        }
    }

    void xtoq(MatrixView& q, const MatrixView& x)
    {
        for(s32 i=0; i<x.rows(); ++i){
            for(s32 j=0; j<x.cols(); ++j){
                q(i,j) /= x(i,i);
            }
            for(s32 j=i+1; j<x.rows(); ++j){
                for(s32 k=0; k<x.cols(); ++k){
                    q(j,k) -= x(j,i) * q(i,k);
                }
            }
        }
    }

    s32 eigen(s32 n, MatrixView& m, VectorView& d, VectorView& e, f32 epsilon, s32 maxIteration)
    {
        LASSERT(0<=n);

        tridiagonalize(n, m, &d[0], &e[1]);
        e[0] = 0.0;
        for(s32 h=n-1; 0<h; --h){
            s32 j=h;
            while(0<j && epsilon * (absolute(d[j-1])+absolute(d[j])) < absolute(e[j])){
                --j;
            }
            if(j==h){
                continue;
            }
            s32 iteration=0;
            element_type e0, e1;
            do{
                if(maxIteration<++iteration){
                    return -1;
                }
                element_type w = (d[h-1] - d[h])*0.5;
                element_type t = e[h]*e[h];
                element_type s = lmath::sqrt(w*w+t);
                if(w<0.0){
                    s = -s;
                }
                element_type x = d[j] - d[h] + t/(w+s);
                element_type y = e[j+1];
                for(s32 k=j; k<h; ++k){
                    element_type c;
                    element_type s;
                    if(absolute(y)<=absolute(x)){
                        t = -y/x;
                        c = 1.0/lmath::sqrt(t*t+1.0);
                        s = t*c;
                    }else{
                        t = -x/y;
                        s = 1.0/lmath::sqrt(t*t+1.0);
                        c = t*s;
                    }
                    w = d[k] - d[k+1];
                    t = (w*s + 2.0*c*e[k+1])*s;
                    d[k] -= t;
                    d[k+1] += t;
                    if(j<k){
                        e[k] = c*e[k] - s*y;
                    }
                    e[k+1] += s*(c*w - 2.0*s*e[k+1]);

                    //固有ベクトルを計算
                    for(s32 i=0; i<n; ++i){
                        x = m(i,k);
                        y = m(i,k+1);
                        m(i,k) = c*x - s*y;
                        m(i,k+1) = s*x + c*y;
                    }
                    if(k<(h-1)){
                        x = e[k+1];
                        y = -s*e[k+2];
                        e[k+2] *= c;
                    }
                }//for(s32 k=j;

                e0 = epsilon*(absolute(d[h-1])+absolute(d[h]));
                e1 = absolute(e[h]);
            } while(e0<=e1);
        }

        //固有値, 固有ベクトルを降順に整列
        for(s32 k=0; k<n-1; ++k){
            s32 h=k;
            element_type t = d[h];
            for(s32 i=k+1; i<n; ++i){
                if(t<d[i]){
                    h=i;
                    t = d[h];
                }
                d[h] = d[k];
                d[k] = t;
                for(s32 j=0; j<m.cols(); ++j){
                    lmath::swap(m(j,h), m(j,k));
                }
            }
        }
        return 0;
    }
#endif
}
