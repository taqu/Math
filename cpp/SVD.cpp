/**
@file SVD.cpp
@author t-sakai
@date 2017/06/02 create
*/
#include "SVD.h"
#include <stdio.h>
#include <limits>
#include "Householder.h"
#include "QR.h"

namespace lmath
{
    DiagonalVector::DiagonalVector()
        :rows_(0)
        ,cols_(0)
        ,x_(NULL)
    {}

    DiagonalVector::DiagonalVector(s32 rows, s32 cols, element_type* x)
        :rows_(rows)
        ,cols_(cols)
        ,x_(x)
    {}

    DiagonalVector::~DiagonalVector()
    {}

    element_type DiagonalVector::operator[](s32 index) const
    {
        LASSERT(0<=index && index<rows_);
        return x_[cols_*index + index];
    }

    element_type& DiagonalVector::operator[](s32 index)
    {
        LASSERT(0<=index && index<rows_);
        return x_[cols_*index + index];
    }


    void bidiagonalization(MatrixView& U, MatrixView& B, MatrixView& V, MatrixView& A)
    {
        LASSERT(B.rows() == A.rows());
        LASSERT(B.cols() == A.cols());

        LASSERT(U.rows() == U.cols());
        LASSERT(V.rows() == V.cols());
        LASSERT(U.rows() == A.rows());
        LASSERT(V.cols() == A.cols());

        s32 size = (A.cols()<A.rows())? A.rows() : A.cols();
        B = A;
        element_type* hbuffer = LNEW element_type[B.rows()*B.rows()];
        element_type* hpbuffer = LNEW element_type[B.cols()*B.cols()];
        element_type* tbuffer = LNEW element_type[size*size];
        element_type* tv0buffer = LNEW element_type[B.rows()];
        element_type* tv1buffer = LNEW element_type[B.cols()];

        U.identity();
        V.identity();
        MatrixView H(B.rows(), B.rows(), hbuffer);
        MatrixView Hprime(B.cols(), B.cols(), hpbuffer);
        MatrixView T0(B.rows(), B.cols(), tbuffer);
        MatrixView T1(B.rows(), B.rows(), tbuffer);
        MatrixView T2(B.cols(), B.cols(), tbuffer);
        VectorView tv0(B.rows(), tv0buffer);
        VectorView tv1(B.cols(), tv1buffer);

        for(s32 k=0; k<B.cols(); ++k){
            //eliminate non-zeros below the diagonal
            tv0.copy(B.getCol(k));
            householder_matrix(H, tv0, k);
            B = mul(T0, H, B);
            U = mul(T1, U, H);
            //eliminate no-zeros to the right to the superdiagonal by working with transpose
            if(k<B.cols()-2){
                tv1.copy(B.getRow(k));
                householder_matrix(Hprime, tv1, k+1);
                transposeSquare(Hprime);
                B = mul(T0, B, Hprime);
                V = mul(T2, Hprime, V);
            }
        }
        LDELETE_ARRAY(tv1buffer);
        LDELETE_ARRAY(tv0buffer);
        LDELETE_ARRAY(tbuffer);
        LDELETE_ARRAY(hpbuffer);
        LDELETE_ARRAY(hbuffer);
    }

    void rotGivens(element_type& c, element_type& s, element_type& r, element_type f, element_type g)
    {
        element_type af = lmath::absolute(f);
        if(af<LMATH_EPSILON){
            c = 0.0;
            s = 1.0;
            r = g;
        }else if(lmath::absolute(g)<af){
            element_type t = g/f;
            element_type t1 = lmath::sqrt(1.0+t*t);
            c = 1.0/t1;
            s = t*c;
            r = f*t1;
        }else{
            element_type t = f/g;
            element_type t1 = lmath::sqrt(1.0+t*t);
            s = 1.0/t1;
            c = t*s;
            r = g*t1;
        }
    }

    void rotGivens(element_type& c, element_type& s, element_type f, element_type g)
    {
        element_type af = lmath::absolute(f);
        if(af<LMATH_EPSILON){
            c = 0.0;
            s = 1.0;
        }else if(lmath::absolute(g)<af){
            element_type t = g/f;
            element_type t1 = lmath::sqrt(1.0+t*t);
            c = 1.0/t1;
            s = t*c;
        }else{
            element_type t = f/g;
            element_type t1 = lmath::sqrt(1.0+t*t);
            s = 1.0/t1;
            c = t*s;
        }
    }

    void msweep(MatrixView& B)
    {
        LASSERT(B.rows() == B.cols());
        element_type* qbuffer = LNEW element_type[B.rows()*B.cols()];
        element_type* tbuffer = LNEW element_type[B.rows()*B.cols()];
        MatrixView Q(B.rows(), B.cols(), qbuffer);
        MatrixView T(B.rows(), B.cols(), tbuffer);

        Q.identity();
        s32 r = B.rows();
        for(s32 i=0; i<r-1; ++i){
            element_type c,s,r;
            rotGivens(c,s,r, B(i,i), B(i,i+1));
            //construct matrix Q and multiply on the right by Q'
            //this annihilates both B(i-1,i+1) and B(i,i+1)
            //but makes B(i+1,i) non-zero
            Q(i,i) = c; Q(i,i+1) = s;
            Q(i+1,i) = -s; Q(i+1, i+1) = c;
            transposeSquare(Q);
            B=mul(T,B,Q);
            Q(i  ,i) = 1.0; Q(i,   i+1) = 0.0;
            Q(i+1,i) = 0.0; Q(i+1, i+1) = 1.0;

            rotGivens(c,s,r, B(i,i), B(i+1,i));
            //construct matrix Q and multiply on the left by Q
            //this annihilates B(i+1,i), but makes B(i,i+1) and
            //B(i,i+2) non-zero
            Q(i,i) = c; Q(i,i+1) = s;
            Q(i+1,i) = -s; Q(i+1, i+1) = c;
            B=mul(T,Q,B);
            Q(i  ,i) = 1.0; Q(i,   i+1) = 0.0;
            Q(i+1,i) = 0.0; Q(i+1, i+1) = 1.0;
        }
        LDELETE_ARRAY(tbuffer);
        LDELETE_ARRAY(qbuffer);
    }

    void msweep2(MatrixView& U, MatrixView& B, MatrixView& V)
    {
        LASSERT(B.rows() == B.cols());
        element_type* qbuffer = LNEW element_type[B.rows()*B.cols()];
        element_type* tbuffer = LNEW element_type[B.rows()*B.cols()];
        MatrixView Q(B.rows(), B.cols(), qbuffer);
        MatrixView T(B.rows(), B.cols(), tbuffer);

        Q.identity();
        s32 r = B.rows();
        for(s32 i=0; i<r-1; ++i){
            element_type c,s,r;
            rotGivens(c,s,r, B(i,i), B(i,i+1));
            //construct matrix Q and multiply on the right by Q'
            //this annihilates both B(i-1,i+1) and B(i,i+1)
            //but makes B(i+1,i) non-zero
            Q(i,i) = c; Q(i,i+1) = -s;
            Q(i+1,i) = s; Q(i+1, i+1) = c;
            B=mul(T,B,Q);
            Q(i,i+1) = s; Q(i+1,i) = -s;
            V=mul(T,Q,V);
            Q(i  ,i) = 1.0; Q(i,   i+1) = 0.0;
            Q(i+1,i) = 0.0; Q(i+1, i+1) = 1.0;

            rotGivens(c,s,r, B(i,i), B(i+1,i));
            //construct matrix Q and multiply on the left by Q
            //this annihilates B(i+1,i), but makes B(i,i+1) and
            //B(i,i+2) non-zero
            Q(i,i) = c; Q(i,i+1) = s;
            Q(i+1,i) = -s; Q(i+1, i+1) = c;
            B=mul(T,Q,B);
            Q(i,i+1) = -s; Q(i+1,i) = s;
            U=mul(T,U,Q);
            Q(i  ,i) = 1.0; Q(i,   i+1) = 0.0;
            Q(i+1,i) = 0.0; Q(i+1, i+1) = 1.0;
        }
        LDELETE_ARRAY(tbuffer);
        LDELETE_ARRAY(qbuffer);
    }

    void msweep2(s32 lower, s32 upper, MatrixView& U, MatrixView& B, MatrixView& V)
    {
        s32 size = (B.rows()<B.cols())? B.cols() : B.rows();
        element_type* qbuffer = LNEW element_type[B.rows()*B.rows() + B.cols()*B.cols()];
        element_type* tbuffer = LNEW element_type[size*size];
        MatrixView Q0(B.rows(), B.rows(), qbuffer);
        MatrixView Q1(B.cols(), B.cols(), qbuffer+B.rows()*B.rows());
        MatrixView T0(B.rows(), B.cols(), tbuffer);
        MatrixView T1(B.rows(), B.rows(), tbuffer);
        MatrixView T2(B.cols(), B.cols(), tbuffer);

        Q0.identity();
        Q1.identity();
        for(s32 i=lower; i<upper; ++i){
            element_type c,s,r;
            rotGivens(c,s,r, B(i,i), B(i,i+1));
            //construct matrix Q and multiply on the right by Q'
            //this annihilates both B(i-1,i+1) and B(i,i+1)
            //but makes B(i+1,i) non-zero
            Q1(i,i) = c; Q1(i,i+1) = -s;
            Q1(i+1,i) = s; Q1(i+1, i+1) = c;
            B=mul(T0,B,Q1);
            Q1(i,i+1) = s; Q1(i+1,i) = -s;
            V=mul(T2,Q1,V);
            Q1(i  ,i) = 1.0; Q1(i,   i+1) = 0.0;
            Q1(i+1,i) = 0.0; Q1(i+1, i+1) = 1.0;

            rotGivens(c,s,r, B(i,i), B(i+1,i));
            //construct matrix Q and multiply on the left by Q
            //this annihilates B(i+1,i), but makes B(i,i+1) and
            //B(i,i+2) non-zero
            Q0(i,i) = c; Q0(i,i+1) = s;
            Q0(i+1,i) = -s; Q0(i+1, i+1) = c;
            B=mul(T0,Q0,B);
            Q0(i,i+1) = -s; Q0(i+1,i) = s;
            U=mul(T1,U,Q0);
            Q0(i  ,i) = 1.0; Q0(i,   i+1) = 0.0;
            Q0(i+1,i) = 0.0; Q0(i+1, i+1) = 1.0;
        }
        LDELETE_ARRAY(tbuffer);
        LDELETE_ARRAY(qbuffer);
    }

    void vsweep(VectorView& d, VectorView& e)
    {
        LASSERT(d.size() ==  (e.size()+1));
        element_type c=1.0, cold=1.0;
        element_type s, sold, r;
        for(s32 i=0; i<e.size(); ++i){
            rotGivens(c,s,r,c*d[i],e[i]);
            if(0!=i){
                e[i-1] = r*sold;
            }
            rotGivens(cold,sold,d[i], cold*r, d[i+1]*s);
        }
        element_type h = c*d[e.size()];
        e[e.size()-1] = h*sold;
        d[d.size()-1] = h*cold;
    }

    void vsweep(s32 lower, s32 upper, VectorView& d, VectorView& e)
    {
        LASSERT(d.size() ==  (e.size()+1));
        LASSERT(upper<=e.size());

        element_type c=1.0, cold=1.0;
        element_type s, sold, r;
        for(s32 i=lower; i<upper; ++i){
            rotGivens(c,s,r,c*d[i],e[i]);
            if(lower!=i){
                e[i-1] = r*sold;
            }
            rotGivens(cold,sold,d[i], cold*r, d[i+1]*s);
        }
        element_type h = c*d[upper];
        e[upper-1] = h*sold;
        d[upper] = h*cold;
    }

    s32 svd(VectorView& d, VectorView& e, element_type epsilon, s32 maxIterationFactor)
    {
        LASSERT(d.size() == (e.size()+1));
        element_type TOL = 100.0*epsilon;
        s32 maxIteration = maxIterationFactor*d.size()*d.size();

        element_type* lambdabuffer = LNEW element_type[d.size()];
        element_type* mubuffer = LNEW element_type[d.size()];

        VectorView lambda(d.size(), lambdabuffer);
        VectorView mu(d.size(), mubuffer);

        //The following convergence criterion is discussed by
        //Demmel and Kahan.  First, estimate the smallest singular value.
        lambda[d.size()-1] = absolute(d[d.size()-1]);
        for(s32 i=d.size()-2; 0<=i; --i){
            lambda[i] = absolute(d[i]) * lambda[i+1]/(lambda[i+1] + absolute(e[i]));
        }

        mu[0] = absolute(d[0]);
        for(s32 i=1;i<d.size(); ++i){
            mu[i] = absolute(d[i]) * mu[i-1]/(mu[i-1]+absolute(e[i-1]));
        }
        element_type sigmaLower = minimum(lambda.minimum(), mu.minimum());
        element_type threshold = maximum(TOL*sigmaLower, maxIteration*std::numeric_limits<element_type>::min());
        LDELETE_ARRAY(mubuffer);
        LDELETE_ARRAY(lambdabuffer);

        s32 iterations;
        s32 iUpper = d.size()-2;
        s32 iLower = 0;
        for(iterations=0; iterations<maxIteration; ++iterations){
            //reduce problem size when some zeros are on the superdiagonal

            //how many zeros are near the bottom right?
            for(s32 i=iUpper; 0<=i; --i){
                iUpper=i;
                if(threshold<absolute(e[i])){
                    break;
                }
            }
            //how many zeros are near the top left?
            {
                s32 j=iUpper;
                for(s32 i=iLower; i<iUpper; ++i){
                    if(threshold<absolute(e[i])){
                        j=i;
                        break;
                    }
                }
                iLower = j;
            }

            if((iUpper == iLower && absolute(e[iUpper])<=threshold) || (iUpper<iLower)){
                //all done, sort singular values in ascending order
                return iterations;
            }

            //do a sweep
            vsweep(iLower, iUpper+1, d, e);
        }
        return -1;
    }

    s32 svd(MatrixView& U, MatrixView& A, MatrixView& V, element_type epsilon, s32 maxIterationFactor)
    {
        //Bidiagonalize
        element_type* bbuffer = LNEW element_type[A.rows()*A.cols()];

        MatrixView B(A.rows(), A.cols(), bbuffer);
        bidiagonalization(U,B,V,A);

        s32 diag = (A.rows()<A.cols())? A.rows() : A.cols();
        DiagonalVector d(diag, A.cols(), &B(0,0));
        DiagonalVector e(diag-1, A.cols(), &B(0,1));

        element_type TOL = 100.0*epsilon;
        s32 maxIteration = maxIterationFactor*d.size()*d.size();

        element_type* lambdabuffer = LNEW element_type[d.size()];
        element_type* mubuffer = LNEW element_type[d.size()];

        VectorView lambda(diag, lambdabuffer);
        VectorView mu(diag, mubuffer);

        //The following convergence criterion is discussed by
        //Demmel and Kahan.  First, estimate the smallest singular value.
        lambda[d.size()-1] = absolute(d[d.size()-1]);
        for(s32 i=d.size()-2; 0<=i; --i){
            lambda[i] = absolute(d[i]) * lambda[i+1]/(lambda[i+1] + absolute(e[i]));
        }

        mu[0] = absolute(d[0]);
        for(s32 i=1;i<d.size(); ++i){
            mu[i] = absolute(d[i]) * mu[i-1]/(mu[i-1]+absolute(e[i-1]));
        }
        element_type sigmaLower = minimum(lambda.minimum(), mu.minimum());
        element_type threshold = maximum(TOL*sigmaLower, maxIteration*std::numeric_limits<element_type>::min());
        LDELETE_ARRAY(mubuffer);
        LDELETE_ARRAY(lambdabuffer);

        s32 iterations;
        s32 iUpper = d.size()-2;
        s32 iLower = 0;
        for(iterations=0; iterations<maxIteration; ++iterations){
            //reduce problem size when some zeros are on the superdiagonal

            //how many zeros are near the bottom right?
            for(s32 i=iUpper; 0<=i; --i){
                iUpper=i;
                if(threshold<absolute(e[i])){
                    break;
                }
            }
            //how many zeros are near the top left?
            {
                s32 j=iUpper;
                for(s32 i=iLower; i<iUpper; ++i){
                    if(threshold<absolute(e[i])){
                        j=i;
                        break;
                    }
                }
                iLower = j;
            }

            if((iUpper == iLower && absolute(e[iUpper])<=threshold) || (iUpper<iLower)){
#if 0
                //all done, sort singular values in ascending order
                for(s32 i=0; i<d.size(); ++i){
                    s32 k=i;
                    for(s32 j=i+1; j<d.size(); ++j){
                        if(lmath::absolute(d[j])<lmath::absolute(d[k])){
                            k=j;
                        }
                    }
                    if(k==i){
                        continue;
                    }
                    lmath::swap(d[k],d[i]);
                    U.swapCols(i,k);
                    V.swapRows(i,k);
                }
#endif
                A = B;
                LDELETE_ARRAY(bbuffer);
                return iterations;
            }

            //do a sweep
            //vsweep(iLower, iUpper+1, d, e);
            msweep2(iLower, iUpper+1, U, B, V);
        }
        LDELETE_ARRAY(bbuffer);
        return -1;
    }

    void pseudoInverse(MatrixView& Aplus, MatrixView& U, MatrixView& A, MatrixView& V)
    {
        LASSERT(Aplus.rows() == A.cols());
        LASSERT(Aplus.cols() == A.rows());
        element_type* tbuffer = LNEW element_type[A.rows()*A.cols()];
        MatrixView T(Aplus.rows(), Aplus.cols(), tbuffer);

        Aplus.zero();
        s32 diag = (A.rows()<A.cols())? A.rows() : A.cols();
        for(s32 i=0; i<diag; ++i){
             element_type a = A(i,i);
             if(LMATH_EPSILON<lmath::absolute(a)){
                 Aplus(i,i) = 1.0/a;
             }
        }

        transposeSquare(V);
        transposeSquare(U);
        Aplus = mul(T,V,Aplus);
        Aplus = mul(T,Aplus,U);
        LDELETE_ARRAY(tbuffer);
    }

#if 0
    s32 svd(MatrixView& u, VectorView& sigma, MatrixView& vt, MatrixView& m, f32 epsilon, s32 maxIteration)
    {
        Matrix tmp_mt(m.rows(), m.cols()); 
        MatrixView mt(tmp_mt);
        m.transpose(mt);

        Vector tmp_sigma1(sigma.size());
        VectorView sigma1(tmp_sigma1);

        Vector tmp_e(sigma.size());
        VectorView e(tmp_e);

        Matrix tmp_ut(u.rows(), u.cols()); 
        MatrixView ut(tmp_ut);

        mul(ut, m, mt);
        s32 ret=0;
        if(0 != eigen(ut.cols(), ut, sigma, e, epsilon, maxIteration)){
            ret |= 0x01U;
        }
        ut.transpose(u);

        mul(vt, mt, m);
        if(0 != eigen(vt.cols(), vt, sigma1, e, epsilon, maxIteration)){
            ret |= 0x02U;
        }
        return ret;
    }

    //-----------------------------------------
    //---
    //--- SVD
    //---
    //-----------------------------------------
    SVD::SVD()
        :cols_(0)
        ,rows_(0)
        ,buffer_(NULL)
    {
    }

    SVD::SVD(s32 cols, s32 rows)
        :cols_(0)
        ,rows_(0)
        ,buffer_(NULL)
    {
        reset(cols, rows);
    }

    SVD::SVD(SVD&& rhs)
        :cols_(rhs.cols_)
        ,rows_(rhs.rows_)
        ,buffer_(rhs.buffer_)
    {
        rhs.cols_ = 0;
        rhs.rows_ = 0;
        rhs.buffer_ = NULL;
    }

    SVD::~SVD()
    {
        cols_ = rows_ = 0;
        LDELETE_ARRAY(buffer_);
    }

    void SVD::reset(s32 cols, s32 rows)
    {
        LASSERT(0<=cols);
        LASSERT(0<=rows);
        if(cols_ == cols && rows_ == rows){
            return;
        }
        s32 size0 = cols*rows;
        s32 size1 = rows*rows;
        s32 size2 = cols*cols;

        s32 n = (rows<cols)? cols : rows;
        if(size0 != (cols_*rows_)){
            LDELETE_ARRAY(buffer_);
            s32 total = size0*3 + size1*2 + size2 + n*2 + cols + rows;
            buffer_ = LNEW element_type[total];
        }
        cols_ = cols;
        rows_ = rows;

        s32 offset=0;

        a_ = MatrixView(cols_, rows_, buffer_+offset);
        offset += size0;

        at_ = MatrixView(rows_, cols_, buffer_+offset);
        offset += size0;

        ap_ = MatrixView(rows_, rows_, buffer_+offset);
        offset += size1;

        ut_ = MatrixView(rows_, rows_, buffer_+offset);
        offset += size1;

        v_ = MatrixView(cols_, cols_, buffer_+offset);
        offset += size2;

        sp_ = MatrixView(rows_, cols_, buffer_+offset);
        offset += size0;

        sigma_ = VectorView(n, buffer_+offset);
        offset += n;

        e_ = VectorView(n, buffer_+offset);
        offset += n;

        x_ = VectorView(cols_, buffer_+offset);
        offset += cols_;

        b_ = VectorView(rows_, buffer_+offset);
        offset += rows_;
    }

    s32 SVD::solve(element_type epsilon, s32 maxIteration, element_type truncate)
    {
        template<typename MatrixType, int QRPreconditioner>
        void JacobiSVD<MatrixType, QRPreconditioner>::allocate(Index rows, Index cols, unsigned int computationOptions)
        {
            eigen_assert(rows >= 0 && cols >= 0);

            if (m_isAllocated &&
                rows == m_rows &&
                cols == m_cols &&
                computationOptions == m_computationOptions)
            {
                return;
            }

            m_rows = rows;
            m_cols = cols;
            m_isInitialized = false;
            m_isAllocated = true;
            m_computationOptions = computationOptions;
            m_computeFullU = (computationOptions & ComputeFullU) != 0;
            m_computeThinU = (computationOptions & ComputeThinU) != 0;
            m_computeFullV = (computationOptions & ComputeFullV) != 0;
            m_computeThinV = (computationOptions & ComputeThinV) != 0;
            eigen_assert(!(m_computeFullU && m_computeThinU) && "JacobiSVD: you can't ask for both full and thin U");
            eigen_assert(!(m_computeFullV && m_computeThinV) && "JacobiSVD: you can't ask for both full and thin V");
            eigen_assert(EIGEN_IMPLIES(m_computeThinU || m_computeThinV, MatrixType::ColsAtCompileTime==Dynamic) &&
                "JacobiSVD: thin U and V are only available when your matrix has a dynamic number of columns.");
            if (QRPreconditioner == FullPivHouseholderQRPreconditioner)
            {
                eigen_assert(!(m_computeThinU || m_computeThinV) &&
                    "JacobiSVD: can't compute thin U or thin V with the FullPivHouseholderQR preconditioner. "
                    "Use the ColPivHouseholderQR preconditioner instead.");
            }
            m_diagSize = (std::min)(m_rows, m_cols);
            m_singularValues.resize(m_diagSize);
            if(RowsAtCompileTime==Dynamic)
                m_matrixU.resize(m_rows, m_computeFullU ? m_rows
                    : m_computeThinU ? m_diagSize
                    : 0);
            if(ColsAtCompileTime==Dynamic)
                m_matrixV.resize(m_cols, m_computeFullV ? m_cols
                    : m_computeThinV ? m_diagSize
                    : 0);
            m_workMatrix.resize(m_diagSize, m_diagSize);

            if(m_cols>m_rows)   m_qr_precond_morecols.allocate(*this);
            if(m_rows>m_cols)   m_qr_precond_morerows.allocate(*this);
            if(m_rows!=m_cols)  m_scaledMatrix.resize(rows,cols);
        }

        template<typename MatrixType, int QRPreconditioner>
        JacobiSVD<MatrixType, QRPreconditioner>&
            JacobiSVD<MatrixType, QRPreconditioner>::compute(const MatrixType& matrix, unsigned int computationOptions)
        {
            using std::abs;
            allocate(matrix.rows(), matrix.cols(), computationOptions);

            // currently we stop when we reach precision 2*epsilon as the last bit of precision can require an unreasonable number of iterations,
            // only worsening the precision of U and V as we accumulate more rotations
            const RealScalar precision = RealScalar(2) * NumTraits<Scalar>::epsilon();

            // limit for denormal numbers to be considered zero in order to avoid infinite loops (see bug 286)
            const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();

            // Scaling factor to reduce over/under-flows
            RealScalar scale = matrix.cwiseAbs().maxCoeff();
            if(scale==RealScalar(0)) scale = RealScalar(1);

            /*** step 1. The R-SVD step: we use a QR decomposition to reduce to the case of a square matrix */

            if(m_rows!=m_cols)
            {
                m_scaledMatrix = matrix / scale;
                m_qr_precond_morecols.run(*this, m_scaledMatrix);
                m_qr_precond_morerows.run(*this, m_scaledMatrix);
            }
            else
            {
                m_workMatrix = matrix.block(0,0,m_diagSize,m_diagSize) / scale;
                if(m_computeFullU) m_matrixU.setIdentity(m_rows,m_rows);
                if(m_computeThinU) m_matrixU.setIdentity(m_rows,m_diagSize);
                if(m_computeFullV) m_matrixV.setIdentity(m_cols,m_cols);
                if(m_computeThinV) m_matrixV.setIdentity(m_cols, m_diagSize);
            }

            /*** step 2. The main Jacobi SVD iteration. ***/
            RealScalar maxDiagEntry = m_workMatrix.cwiseAbs().diagonal().maxCoeff();

            bool finished = false;
            while(!finished)
            {
                finished = true;

                // do a sweep: for all index pairs (p,q), perform SVD of the corresponding 2x2 sub-matrix

                for(Index p = 1; p < m_diagSize; ++p)
                {
                    for(Index q = 0; q < p; ++q)
                    {
                        // if this 2x2 sub-matrix is not diagonal already...
                        // notice that this comparison will evaluate to false if any NaN is involved, ensuring that NaN's don't
                        // keep us iterating forever. Similarly, small denormal numbers are considered zero.
                        RealScalar threshold = numext::maxi<RealScalar>(considerAsZero, precision * maxDiagEntry);
                        if(abs(m_workMatrix.coeff(p,q))>threshold || abs(m_workMatrix.coeff(q,p)) > threshold)
                        {
                            finished = false;
                            // perform SVD decomposition of 2x2 sub-matrix corresponding to indices p,q to make it diagonal
                            // the complex to real operation returns true if the updated 2x2 block is not already diagonal
                            if(internal::svd_precondition_2x2_block_to_be_real<MatrixType, QRPreconditioner>::run(m_workMatrix, *this, p, q, maxDiagEntry))
                            {
                                JacobiRotation<RealScalar> j_left, j_right;
                                internal::real_2x2_jacobi_svd(m_workMatrix, p, q, &j_left, &j_right);

                                // accumulate resulting Jacobi rotations
                                m_workMatrix.applyOnTheLeft(p,q,j_left);
                                if(computeU()) m_matrixU.applyOnTheRight(p,q,j_left.transpose());

                                m_workMatrix.applyOnTheRight(p,q,j_right);
                                if(computeV()) m_matrixV.applyOnTheRight(p,q,j_right);

                                // keep track of the largest diagonal coefficient
                                maxDiagEntry = numext::maxi<RealScalar>(maxDiagEntry,numext::maxi<RealScalar>(abs(m_workMatrix.coeff(p,p)), abs(m_workMatrix.coeff(q,q))));
                            }
                        }
                    }
                }
            }

            /*** step 3. The work matrix is now diagonal, so ensure it's positive so its diagonal entries are the singular values ***/

            for(Index i = 0; i < m_diagSize; ++i)
            {
                // For a complex matrix, some diagonal coefficients might note have been
                // treated by svd_precondition_2x2_block_to_be_real, and the imaginary part
                // of some diagonal entry might not be null.
                if(NumTraits<Scalar>::IsComplex && abs(numext::imag(m_workMatrix.coeff(i,i)))>considerAsZero)
                {
                    RealScalar a = abs(m_workMatrix.coeff(i,i));
                    m_singularValues.coeffRef(i) = abs(a);
                    if(computeU()) m_matrixU.col(i) *= m_workMatrix.coeff(i,i)/a;
                }
                else
                {
                    // m_workMatrix.coeff(i,i) is already real, no difficulty:
                    RealScalar a = numext::real(m_workMatrix.coeff(i,i));
                    m_singularValues.coeffRef(i) = abs(a);
                    if(computeU() && (a<RealScalar(0))) m_matrixU.col(i) = -m_matrixU.col(i);
                }
            }

            m_singularValues *= scale;

            /*** step 4. Sort singular values in descending order and compute the number of nonzero singular values ***/

            m_nonzeroSingularValues = m_diagSize;
            for(Index i = 0; i < m_diagSize; i++)
            {
                Index pos;
                RealScalar maxRemainingSingularValue = m_singularValues.tail(m_diagSize-i).maxCoeff(&pos);
                if(maxRemainingSingularValue == RealScalar(0))
                {
                    m_nonzeroSingularValues = i;
                    break;
                }
                if(pos)
                {
                    pos += i;
                    std::swap(m_singularValues.coeffRef(i), m_singularValues.coeffRef(pos));
                    if(computeU()) m_matrixU.col(pos).swap(m_matrixU.col(i));
                    if(computeV()) m_matrixV.col(pos).swap(m_matrixV.col(i));
                }
            }

            m_isInitialized = true;
            return *this;
        }




        a_.transpose(at_);

        s32 ret=0;
        mul(ut_, a_, at_);
        if(0 != eigen(ut_.cols(), ut_, sigma_, e_, epsilon, maxIteration)){
            ret |= 0x01U;
        }

        mul(v_, at_, a_);
        if(0 != eigen(v_.cols(), v_, sigma_, e_, epsilon, maxIteration)){
            ret |= 0x02U;
        }
        v_.transpose();

        sp_.setIdentity();
        for(s32 i=0; i<cols_; ++i){
            if(truncate<sigma_[i]){
                sp_(i,i) = 1.0f/lmath::sqrt(sigma_[i]);
            }else{
                sp_(i,i) = 0.0;
            }
        }
        mul(at_, sp_, ut_);
        mul(ap_, v_, at_);
        for(s32 i=0; i<cols_; ++i){
            x_[i] = 0.0;
            for(s32 j=0; j<rows_; ++j){
                x_[i] += ap_(j,i) * b_[j];
            }
        }
        return ret;
    }

    SVD& SVD::operator=(SVD&& rhs)
    {
        cols_ = rhs.cols_;
        rows_ = rhs.rows_;
        buffer_ = rhs.buffer_;

        rhs.cols_ = 0;
        rhs.rows_ = 0;
        rhs.buffer_ = NULL;
        return *this;
    }
#endif
}
