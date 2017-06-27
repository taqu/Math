/**
@file Matrix.cpp
@author t-sakai
@date 2017/06/02 create
*/
#include "Matrix.h"
#include <stdio.h>
#include <memory.h>

namespace lmath
{
    //--------------------------------------------------
    //---
    //--- MatrixView
    //---
    //--------------------------------------------------
    MatrixView::MatrixView()
        :rows_(0)
        ,cols_(0)
        ,x_(NULL)
    {
    }

    MatrixView::MatrixView(s32 rows, s32 cols, element_type* x)
        :rows_(rows)
        ,cols_(cols)
        ,x_(x)
    {
        LASSERT(0<=rows);
        LASSERT(0<=cols);
        LASSERT(NULL != x);
    }

    MatrixView::MatrixView(MatrixView&& m)
        :rows_(m.rows_)
        ,cols_(m.cols_)
        ,x_(m.x_)
    {
        m.rows_ = 0;
        m.cols_ = 0;
        m.x_ = NULL;
    }

    MatrixView::~MatrixView()
    {
        cols_ = 0;
        rows_ = 0;
        x_ = NULL;
    }

    MatrixView& MatrixView::copy(const MatrixView& m)
    {
        LASSERT(rows_ == m.rows());
        LASSERT(cols_ == m.cols());
        s32 size = cols_*rows_;
        memcpy(x_, m.x_, sizeof(element_type)*size);
        return *this;
    }

    MatrixView& MatrixView::operator=(MatrixView&& m)
    {
        rows_ = m.rows_;
        cols_ = m.cols_;
        x_ = m.x_;
        m.rows_ = 0;
        m.cols_ = 0;
        m.x_ = NULL;
        return *this;
    }

    void MatrixView::zero()
    {
        memset(x_, 0, sizeof(element_type)*rows_*cols_);
    }

    void MatrixView::identity()
    {
        memset(x_, 0, sizeof(element_type)*rows_*cols_);
        s32 n = (cols_<rows_)? cols_ : rows_;
        for(s32 i=0; i<n; ++i){
            (*this)(i,i) = 1.0;
        }
    }

    void MatrixView::transpose(MatrixView& dst) const
    {
        LASSERT(rows_ == dst.cols());
        LASSERT(cols_ == dst.rows());

        for(s32 i=0; i<rows_; ++i){
            for(s32 j=0; j<cols_; ++j){
                dst(j,i) = (*this)(i,j);
            }
        }
    }

    VectorView MatrixView::getRow(s32 row)
    {
        s32 offset = row*cols_;
        return VectorView(cols_, x_+offset);
    }

    VectorStepView MatrixView::getCol(s32 col)
    {
        LASSERT(0<=col && col<cols());
        return VectorStepView(rows_, cols_, x_+col);
    }

    void MatrixView::swapRows(s32 r0, s32 r1)
    {
        element_type* x0 = x_ + r0*cols();
        element_type* x1 = x_ + r1*cols();
        for(s32 i=0; i<cols(); ++i){
            lmath::swap(x0[i], x1[i]);
        }
    }

    void MatrixView::swapCols(s32 c0, s32 c1)
    {
        element_type* x0 = x_ + c0;
        element_type* x1 = x_ + c1;
        for(s32 i=0; i<rows(); ++i){
            lmath::swap(*x0, *x1);
            x0 += cols();
            x1 += cols();
        }
    }

    MatrixView& mul(MatrixView& dst, const MatrixView& m0, const MatrixView& m1)
    {
        LASSERT(dst.rows() == m0.rows());
        LASSERT(dst.cols() == m1.cols());
        LASSERT(m0.cols() == m1.rows());
        LASSERT(&dst != &m0);
        LASSERT(&dst != &m1);

        for(s32 i=0; i<m0.rows(); ++i){
            for(s32 j=0; j<m1.cols(); ++j){
                element_type t=0.0;
                for(s32 k=0; k<m0.cols(); ++k){
                    t += m0(i,k) * m1(k,j);
                }
                dst(i,j) = t;
            }
        }
        return dst;
    }

    void print(const MatrixView& m)
    {
        for(s32 i=0; i<m.rows(); ++i){
            for(s32 j=0; j<m.cols(); ++j){
                printf("%f, ", m(i,j));
            }
            printf("\n");
        }
    }

    VectorView& mul(VectorView& dst, const VectorView& v, const MatrixView& m)
    {
        LASSERT(dst.size() == m.cols());
        LASSERT(v.size() == m.rows());
        LASSERT(&dst != &v);
        for(s32 i=0; i<dst.size(); ++i){
            element_type t = 0.0;
            for(s32 j=0; j<v.size(); ++j){
                t += v[j]*m(j,i);
            }
            dst[i] = t;
        }
        return dst;
    }

    VectorView& mul(VectorView& dst, const MatrixView& m, const VectorView& v)
    {
        LASSERT(dst.size() == m.rows());
        LASSERT(m.cols() == v.size());
        LASSERT(&dst != &v);
        for(s32 i=0; i<dst.size(); ++i){
            element_type t = 0.0;
            for(s32 j=0; j<v.size(); ++j){
                t += m(i,j)*v[j];
            }
            dst[i] = t;
        }
        return dst;
    }

    void transposeSquare(MatrixView& m)
    {
        LASSERT(m.cols() == m.rows());
        for(s32 i=1; i<m.rows(); ++i){
            for(s32 j=0; j<i; ++j){
                lmath::swap(m(i,j), m(j,i));
            }
        }
    }

}
