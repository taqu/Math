#ifndef INC_LMATH_MATRIX_H__
#define INC_LMATH_MATRIX_H__
/**
@file Matrix.h
@author t-sakai
@date 2017/06/02 create
*/
#include "lmath.h"
#include "Vector.h"

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
        Matrix(s32 rows, s32 cols);
        Matrix(const Matrix& m);
        Matrix(Matrix&& m);
        ~Matrix();

        inline s32 rows() const;
        inline s32 cols() const;
        inline element_type operator()(s32 row, s32 col) const;
        inline element_type& operator()(s32 row, s32 col);
        Matrix& operator=(const Matrix& m)
        {
            return copy(m);
        }
        Matrix& operator=(Matrix&& m);

        Matrix& copy(const Matrix& m);

        void identity();
        void transpose(Matrix& dst) const;
        VectorView getRow(s32 row);
        VectorStepView getCol(s32 col);
    private:
        friend class MatrixView;
        element_type operator[](s32 index) const = delete;
        element_type& operator[](s32 index) = delete;

        s32 rows_;
        s32 cols_;
        element_type* x_;
    };

    inline s32 Matrix::rows() const
    {
        return rows_;
    }

    inline s32 Matrix::cols() const
    {
        return cols_;
    }

    inline Matrix::element_type Matrix::operator()(s32 row, s32 col) const
    {
        LASSERT(0<=row && row<rows_);
        LASSERT(0<=col && col<cols_);
        return x_[row*cols_+col];
    }

    inline Matrix::element_type& Matrix::operator()(s32 row, s32 col)
    {
        LASSERT(0<=row && row<rows_);
        LASSERT(0<=col && col<cols_);
        return x_[row*cols_+col];
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
        MatrixView(s32 rows, s32 cols, element_type* x);
        MatrixView(MatrixView&& m);
        MatrixView(Matrix& m);
        ~MatrixView();

        inline s32 rows() const;
        inline s32 cols() const;
        inline element_type operator()(s32 row, s32 col) const;
        inline element_type& operator()(s32 row, s32 col);
        MatrixView& operator=(const MatrixView& m)
        {
            return copy(m);
        }
        MatrixView& operator=(MatrixView&& m);
        MatrixView& operator=(Matrix& m);

        MatrixView& copy(const MatrixView& m);

        void identity();
        void transpose(MatrixView& dst) const;
        VectorView getRow(s32 row);
        VectorStepView getCol(s32 col);
    private:
        MatrixView(const MatrixView&) = delete;
        element_type operator[](s32 index) const = delete;
        element_type& operator[](s32 index) = delete;

        s32 rows_;
        s32 cols_;
        element_type* x_;
    };

    inline s32 MatrixView::rows() const
    {
        return rows_;
    }

    inline s32 MatrixView::cols() const
    {
        return cols_;
    }

    inline MatrixView::element_type MatrixView::operator()(s32 row, s32 col) const
    {
        LASSERT(0<=row && row<rows_);
        LASSERT(0<=col && col<cols_);
        return x_[row*cols_+col];
    }

    inline MatrixView::element_type& MatrixView::operator()(s32 row, s32 col)
    {
        LASSERT(0<=row && row<rows_);
        LASSERT(0<=col && col<cols_);
        return x_[row*cols_+col];
    }

    MatrixView& mul(MatrixView& dst, const MatrixView& m0, const MatrixView& m1);
    void print(const MatrixView& m);

    //--------------------------------------------------
    VectorView& mul(VectorView& dst, const VectorView& v, const MatrixView& m);
    VectorView& mul(VectorView& dst, const MatrixView& m, const VectorView& v);
    void transposeSquare(MatrixView& m);
}
#endif //INC_LMATH_MATRIX_H__
