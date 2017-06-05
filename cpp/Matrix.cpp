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
        :cols_(0)
        ,rows_(0)
        ,x_(NULL)
    {
    }

    Matrix::Matrix(s32 cols, s32 rows)
        :cols_(cols)
        ,rows_(rows)
        ,x_(NULL)
    {
        LASSERT(0<=cols_);
        LASSERT(0<=rows_);
        x_ = LNEW element_type[cols_*rows_];
    }

    Matrix::Matrix(const Matrix& m)
        :cols_(m.cols_)
        ,rows_(m.rows_)
        ,x_(NULL)
    {
        s32 size = cols_*rows_;
        x_ = LNEW element_type[size];
        memcpy(x_, m.x_, sizeof(element_type)*size);
    }

    Matrix::Matrix(Matrix&& m)
        :cols_(m.cols_)
        ,rows_(m.rows_)
        ,x_(m.x_)
    {
        m.cols_ = 0;
        m.rows_ = 0;
        m.x_ = NULL;
    }

    Matrix::~Matrix()
    {
        cols_ = 0;
        rows_ = 0;
        LDELETE_ARRAY(x_);
    }

    Matrix& Matrix::operator=(const Matrix& m)
    {
        LDELETE_ARRAY(x_);
        cols_ = m.cols_;
        rows_ = m.rows_;
        s32 size = cols_*rows_;
        x_ = LNEW element_type[size];
        memcpy(x_, m.x_, sizeof(element_type)*size);
        return *this;
    }

    Matrix& Matrix::operator=(Matrix&& m)
    {
        LDELETE_ARRAY(x_);
        cols_ = m.cols_;
        rows_ = m.rows_;
        x_ = m.x_;
        m.cols_ = 0;
        m.rows_ = 0;
        m.x_ = NULL;
        return *this;
    }

    void Matrix::setIdentity()
    {
        memset(x_, 0, sizeof(element_type)*cols_*rows_);
        s32 n = (cols_<rows_)? cols_ : rows_;
        for(s32 i=0; i<n; ++i){
            (*this)(i,i) = 1.0;
        }
    }

    void Matrix::transpose(Matrix& dst) const
    {
        for(s32 i=0; i<rows_; ++i){
            for(s32 j=0; j<cols_; ++j){
                dst(i, j) = (*this)(j, i);
            }
        }
    }

    void mul(Matrix& dst, const Matrix& m0, const Matrix& m1)
    {
        for(s32 i=0; i<m0.rows(); ++i){
            for(s32 j=0; j<m1.cols(); ++j){
                element_type t=0.0;
                for(s32 k=0; k<m0.cols(); ++k){
                    t += m0(k,i) * m1(j,k);
                }
                dst(j,i) = t;
            }
        }
    }

    void print(const Matrix& m)
    {
        for(s32 i=0; i<m.rows(); ++i){
            for(s32 j=0; j<m.cols(); ++j){
                printf("%f, ", m[m.getIndex(j,i)]);
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
        :cols_(0)
        ,rows_(0)
        ,x_(NULL)
    {
    }

    MatrixView::MatrixView(s32 cols, s32 rows, element_type* x)
        :cols_(cols)
        ,rows_(rows)
        ,x_(x)
    {
        LASSERT(0<=cols);
        LASSERT(0<=rows);
        LASSERT(NULL != x);
    }

    MatrixView::MatrixView(MatrixView&& m)
        :cols_(m.cols_)
        ,rows_(m.rows_)
        ,x_(m.x_)
    {
        m.cols_ = 0;
        m.rows_ = 0;
        m.x_ = NULL;
    }

    MatrixView::MatrixView(Matrix& m)
        :cols_(m.cols_)
        ,rows_(m.rows_)
        ,x_(m.x_)
    {
    }

    MatrixView::~MatrixView()
    {
        cols_ = 0;
        rows_ = 0;
    }

    MatrixView& MatrixView::operator=(MatrixView&& m)
    {
        cols_ = m.cols_;
        rows_ = m.rows_;
        x_ = m.x_;
        m.cols_ = 0;
        m.rows_ = 0;
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

    void MatrixView::setSize(s32 cols, s32 rows)
    {
        LASSERT((cols*rows)<=(cols_*rows_));
        cols_ = cols;
        rows_ = rows;
    }

    void MatrixView::setIdentity()
    {
        memset(x_, 0, sizeof(element_type)*cols_*rows_);
        s32 n = (cols_<rows_)? cols_ : rows_;
        for(s32 i=0; i<n; ++i){
            (*this)(i,i) = 1.0;
        }
    }

    void MatrixView::transpose()
    {
        for(s32 i=0; i<rows_; ++i){
            for(s32 j=0; j<cols_; ++j){
                lmath::swap((*this)(i,j), (*this)(j,i));
            }
        }
    }


    void MatrixView::transpose(MatrixView& dst) const
    {
        for(s32 i=0; i<rows_; ++i){
            for(s32 j=0; j<cols_; ++j){
                dst(i,j) = (*this)(j,i);
            }
        }
    }
    void mul(MatrixView& dst, const MatrixView& m0, const MatrixView& m1)
    {
        for(s32 i=0; i<m0.rows(); ++i){
            for(s32 j=0; j<m1.cols(); ++j){
                element_type t=0.0;
                for(s32 k=0; k<m0.cols(); ++k){
                    t += m0(k,i) * m1(j,k);
                }
                dst(j,i) = t;
            }
        }
    }

    void print(const MatrixView& m)
    {
        for(s32 i=0; i<m.rows(); ++i){
            for(s32 j=0; j<m.cols(); ++j){
                printf("%f, ", m[m.getIndex(j,i)]);
            }
            printf("\n");
        }
    }
}
