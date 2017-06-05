#ifndef INC_LMATH_MATRIX_H__
#define INC_LMATH_MATRIX_H__
/**
@file Matrix.h
@author t-sakai
@date 2017/06/02 create
*/
#include "lmath.h"

namespace lmath
{
    //--------------------------------------------------
    //---
    //--- Matrix
    //---
    //--------------------------------------------------
    class Matrix
    {
    public:
        typedef lmath::element_type element_type;

        Matrix();
        Matrix(s32 cols, s32 rows);
        Matrix(const Matrix& m);
        Matrix(Matrix&& m);
        ~Matrix();

        inline s32 cols() const;
        inline s32 rows() const;
        inline s32 getIndex(s32 c, s32 r) const;
        inline element_type operator[](s32 index) const;
        inline element_type& operator[](s32 index);
        inline element_type operator()(s32 c, s32 r) const;
        inline element_type& operator()(s32 c, s32 r);
        Matrix& operator=(const Matrix& m);
        Matrix& operator=(Matrix&& m);

        void setIdentity();
        void transpose(Matrix& dst) const;
    private:
        friend class MatrixView;

        s32 cols_;
        s32 rows_;
        element_type* x_;
    };

    inline s32 Matrix::cols() const
    {
        return cols_;
    }

    inline s32 Matrix::rows() const
    {
        return rows_;
    }

    inline s32 Matrix::getIndex(s32 c, s32 r) const
    {
        LASSERT(0<=c && c<cols_);
        LASSERT(0<=r && r<rows_);
        return r*cols_+c;
    }

    inline Matrix::element_type Matrix::operator[](s32 index) const
    {
        LASSERT(0<=index && index<(rows_*cols_));
        return x_[index];
    }

    inline Matrix::element_type& Matrix::operator[](s32 index)
    {
        LASSERT(0<=index && index<(rows_*cols_));
        return x_[index];
    }

    inline Matrix::element_type Matrix::operator()(s32 c, s32 r) const
    {
        LASSERT(0<=c && c<cols_);
        LASSERT(0<=r && r<rows_);
        return x_[r*cols_+c];
    }

    inline Matrix::element_type& Matrix::operator()(s32 c, s32 r)
    {
        LASSERT(0<=c && c<cols_);
        LASSERT(0<=r && r<rows_);
        return x_[r*cols_+c];
    }

    void mul(Matrix& dst, const Matrix& m0, const Matrix& m1);
    void print(const Matrix& m);

    //--------------------------------------------------
    //---
    //--- MatrixView
    //---
    //--------------------------------------------------
    class MatrixView
    {
    public:
        typedef lmath::element_type element_type;

        MatrixView();
        MatrixView(s32 cols, s32 rows, element_type* x);
        MatrixView(MatrixView&& m);
        MatrixView(Matrix& m);
        ~MatrixView();

        inline s32 cols() const;
        inline s32 rows() const;
        inline s32 getIndex(s32 c, s32 r) const;
        inline element_type operator[](s32 index) const;
        inline element_type& operator[](s32 index);
        inline element_type operator()(s32 c, s32 r) const;
        inline element_type& operator()(s32 c, s32 r);
        MatrixView& operator=(MatrixView&& m);
        MatrixView& operator=(Matrix& m);

        void setSize(s32 cols, s32 rows);
        void setIdentity();
        void transpose();
        void transpose(MatrixView& dst) const;
    private:
        MatrixView(const MatrixView&) = delete;
        MatrixView& operator=(const MatrixView&) = delete;

        s32 cols_;
        s32 rows_;
        element_type* x_;
    };

    inline s32 MatrixView::cols() const
    {
        return cols_;
    }

    inline s32 MatrixView::rows() const
    {
        return rows_;
    }

    inline s32 MatrixView::getIndex(s32 c, s32 r) const
    {
        LASSERT(0<=c && c<cols_);
        LASSERT(0<=r && r<rows_);
        return r*cols_+c;
    }

    inline MatrixView::element_type MatrixView::operator[](s32 index) const
    {
        LASSERT(0<=index && index<(rows_*cols_));
        return x_[index];
    }

    inline MatrixView::element_type& MatrixView::operator[](s32 index)
    {
        LASSERT(0<=index && index<(rows_*cols_));
        return x_[index];
    }

    inline MatrixView::element_type MatrixView::operator()(s32 c, s32 r) const
    {
        LASSERT(0<=c && c<cols_);
        LASSERT(0<=r && r<rows_);
        return x_[r*cols_+c];
    }

    inline MatrixView::element_type& MatrixView::operator()(s32 c, s32 r)
    {
        LASSERT(0<=c && c<cols_);
        LASSERT(0<=r && r<rows_);
        return x_[r*cols_+c];
    }

    void mul(MatrixView& dst, const MatrixView& m0, const MatrixView& m1);
    void print(const MatrixView& m);
}
#endif //INC_LMATH_MATRIX_H__
