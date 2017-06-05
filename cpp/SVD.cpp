/**
@file SVD.cpp
@author t-sakai
@date 2017/06/02 create
*/
#include "SVD.h"
#include <stdio.h>
#include "QR.h"

namespace lmath
{
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
            s32 total = size0*4 + size1 + size2 + n*2 + cols + rows;
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
        offset += size0;

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
}
