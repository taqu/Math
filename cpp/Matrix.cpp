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
    //--- Matrix
    //---
    //--------------------------------------------------
    Matrix::Matrix()
        :rows_(0)
        ,cols_(0)
        ,x_(NULL)
    {
    }

    Matrix::Matrix(s32 rows, s32 cols)
        :rows_(rows)
        ,cols_(cols)
        ,x_(NULL)
    {
        LASSERT(0<=rows_);
        LASSERT(0<=cols_);
        x_ = LNEW element_type[rows_*cols_];
    }

    Matrix::Matrix(const Matrix& m)
        :rows_(m.rows_)
        ,cols_(m.cols_)
        ,x_(NULL)
    {
        s32 size = rows_*cols_;
        x_ = LNEW element_type[size];
        memcpy(x_, m.x_, sizeof(element_type)*size);
    }

    Matrix::Matrix(Matrix&& m)
        :rows_(m.rows_)
        ,cols_(m.cols_)
        ,x_(m.x_)
    {
        m.rows_ = 0;
        m.cols_ = 0;
        m.x_ = NULL;
    }

    Matrix::~Matrix()
    {
        rows_ = 0;
        cols_ = 0;
        LDELETE_ARRAY(x_);
    }

    Matrix& Matrix::copy(const Matrix& m)
    {
        LASSERT(m.rows() == rows_);
        LASSERT(m.cols() == cols_);
        s32 size = cols_*rows_;
        memcpy(x_, m.x_, sizeof(element_type)*size);
        return *this;
    }

    Matrix& Matrix::operator=(Matrix&& m)
    {
        LDELETE_ARRAY(x_);
        rows_ = m.rows_;
        cols_ = m.cols_;
        x_ = m.x_;
        m.cols_ = 0;
        m.rows_ = 0;
        m.x_ = NULL;
        return *this;
    }

    void Matrix::identity()
    {
        memset(x_, 0, sizeof(element_type)*cols_*rows_);
        s32 n = (cols_<rows_)? cols_ : rows_;
        for(s32 i=0; i<n; ++i){
            (*this)(i,i) = 1.0;
        }
    }

    void Matrix::transpose(Matrix& dst) const
    {
        LASSERT(dst.rows() == cols_);
        LASSERT(dst.cols() == rows_);
        for(s32 i=0; i<rows_; ++i){
            for(s32 j=0; j<cols_; ++j){
                dst(j, i) = (*this)(i, j);
            }
        }
    }

    VectorView Matrix::getRow(s32 row)
    {
        LASSERT(0<=row && row<rows());
        s32 offset = row*cols_;
        return VectorView(cols_, x_+offset);
    }

    VectorStepView Matrix::getCol(s32 col)
    {
        LASSERT(0<=col && col<cols());
        return VectorStepView(rows_, cols_, x_+col);
    }

    void mul(Matrix& dst, const Matrix& m0, const Matrix& m1)
    {
        LASSERT(m0.cols() == m1.rows());
        LASSERT(dst.rows() == m0.rows());
        LASSERT(dst.cols() == m1.cols());
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
    }

    void print(const Matrix& m)
    {
        for(s32 i=0; i<m.rows(); ++i){
            for(s32 j=0; j<m.cols(); ++j){
                printf("%f, ", m(i,j));
            }
            printf("\n");
        }
    }

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

    MatrixView::MatrixView(Matrix& m)
        :rows_(m.rows_)
        ,cols_(m.cols_)
        ,x_(m.x_)
    {
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

    MatrixView& MatrixView::operator=(Matrix& m)
    {
        cols_ = m.cols_;
        rows_ = m.rows_;
        x_ = m.x_;
        return *this;
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
