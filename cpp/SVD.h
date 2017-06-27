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
    class DiagonalVector
    {
    public:
        DiagonalVector();
        DiagonalVector(s32 rows, s32 cols, element_type* x);
        ~DiagonalVector();

        inline s32 size() const{return rows_;}
        element_type operator[](s32 index) const;
        element_type& operator[](s32 index);
    private:
        DiagonalVector(const DiagonalVector&) = delete;
        DiagonalVector& operator=(const DiagonalVector&) = delete;

        s32 rows_;
        s32 cols_;
        element_type* x_;
    };

    void bidiagonalization(MatrixView& U, MatrixView& B, MatrixView& V, MatrixView& A);

    void rotGivens(element_type& c, element_type& s, element_type& r, element_type f, element_type g);
    void rotGivens(element_type& c, element_type& s, element_type f, element_type g);

    void msweep(MatrixView& B);
    void msweep2(MatrixView& U, MatrixView& B, MatrixView& V);
    void msweep2(s32 lower, s32 upper, MatrixView& U, MatrixView& B, MatrixView& V);
    void vsweep(VectorView& d, VectorView& e);
    void vsweep(s32 lower, s32 upper, VectorView& d, VectorView& e);

    /**
    @return ê¨å˜Ç»ÇÁ0à»è„
    */
    s32 svd(VectorView& d, VectorView& e, element_type epsilon=LMATH_EPSILON, s32 maxIterationFactor=500);

    /**
    @return ê¨å˜Ç»ÇÁà»è„
    */
    s32 svd(MatrixView& U, MatrixView& A, MatrixView& V, element_type epsilon=LMATH_EPSILON, s32 maxIterationFactor=500);

    void pseudoInverse(MatrixView& Aplus, MatrixView& U, MatrixView& A, MatrixView& V);

 }
#endif //INC_LMATH_SVD_H__
