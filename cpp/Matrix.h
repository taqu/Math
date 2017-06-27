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

        MatrixView& copy(const MatrixView& m);

        void zero();
        void identity();
        void transpose(MatrixView& dst) const;
        VectorView getRow(s32 row);
        VectorStepView getCol(s32 col);

        void swapRows(s32 r0, s32 r1);
        void swapCols(s32 c0, s32 c1);
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
