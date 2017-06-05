#ifndef INC_LMATH_SVD_H__
#define INC_LMATH_SVD_H__
/**
@file SVD.h
@author t-sakai
@date 2017/06/02 create
*/
#include "lmath.h"
#include "Vector.h"
#include "Matrix.h"

namespace lmath
{
    /**
    @return 0Ç»ÇÁê¨å˜
    @param u ... r x r
    @param sigma ... max(r,c)
    @param vt ... c x c
    @param m ... r x c
    */
    s32 svd(MatrixView& u, VectorView& sigma, MatrixView& vt, MatrixView& m, f32 epsilon=LMATH_EPSILON, s32 maxIteration=10);

    //-----------------------------------------
    //---
    //--- SVD
    //---
    //-----------------------------------------
    class SVD
    {
    public:
        SVD();
        SVD(s32 cols, s32 rows);
        SVD(SVD&& rhs);
        ~SVD();

        void reset(s32 cols, s32 rows);
        s32 solve(element_type epsilon=LMATH_EPSILON, s32 maxIteration=10, element_type truncate=LMATH_TRUNCATE_EPSILON);

        inline const MatrixView& a() const;
        inline MatrixView& a();
        inline const MatrixView& ut() const;
        inline MatrixView& ut();
        inline const MatrixView& v() const;
        inline MatrixView& v();
        inline const VectorView& sigma() const;
        inline VectorView& sigma();

        inline const VectorView& x() const;
        inline VectorView& x();

        inline const VectorView& b() const;
        inline VectorView& b();

        SVD& operator=(SVD&& rhs);
    private:
        SVD(const SVD&) = delete;
        SVD& operator=(const SVD&) = delete;

        s32 cols_;
        s32 rows_;
        element_type* buffer_;

        MatrixView a_;
        MatrixView at_;
        MatrixView ap_;
        MatrixView ut_;
        MatrixView v_;
        MatrixView sp_;
        VectorView sigma_;
        VectorView e_;

        VectorView x_;
        VectorView b_;
    };

    inline const MatrixView& SVD::a() const
    {
        return a_;
    }

    inline MatrixView& SVD::a()
    {
        return a_;
    }

    inline const MatrixView& SVD::ut() const
    {
        return ut_;
    }

    inline MatrixView& SVD::ut()
    {
        return ut_;
    }

    inline const MatrixView& SVD::v() const
    {
        return v_;
    }

    inline MatrixView& SVD::v()
    {
        return v_;
    }

    inline const VectorView& SVD::sigma() const
    {
        return sigma_;
    }

    inline VectorView& SVD::sigma()
    {
        return sigma_;
    }

    inline const VectorView& SVD::x() const
    {
        return x_;
    }

    inline VectorView& SVD::x()
    {
        return x_;
    }

    inline const VectorView& SVD::b() const
    {
        return b_;
    }

    inline VectorView& SVD::b()
    {
        return b_;
    }
}
#endif //INC_LMATH_SVD_H__
