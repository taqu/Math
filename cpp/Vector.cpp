/**
@file Vector.cpp
@author t-sakai
@date 2017/06/02 create
*/
#include "Vector.h"
#include <memory.h>
#include <stdio.h>

namespace lmath
{
    element_type innerproduct(s32 size, const element_type* v0, const element_type* v1)
    {
        LASSERT(NULL != v0);
        LASSERT(NULL != v1);
        element_type norm=0.0;
        for(s32 i=0; i<size; ++i){
            norm += v0[i]*v1[i];
        }
        return norm;
    }

    //--------------------------------------------------
    //---
    //--- Vector
    //---
    //--------------------------------------------------
    Vector::Vector()
        :size_(0)
        ,x_(NULL)
    {
    }

    Vector::Vector(s32 n)
        :size_(n)
        ,x_(NULL)
    {
        LASSERT(0<=size_);
        x_ = LNEW element_type[size_];
    }

    Vector::Vector(const Vector& v)
        :size_(v.size_)
        ,x_(NULL)
    {
        x_ = LNEW element_type[size_];
        memcpy(x_, v.x_, sizeof(element_type)*size_);
    }

    Vector::Vector(Vector&& v)
        :size_(v.size_)
        ,x_(v.x_)
    {
        v.size_ = 0;
        v.x_ = NULL;
    }

    Vector::~Vector()
    {
        size_ = 0;
        LDELETE_ARRAY(x_);
    }

    Vector& Vector::operator=(const Vector& v)
    {
        LASSERT(size_ == v.size());
        memcpy(x_, v.x_, sizeof(element_type)*size_);
        return *this;
    }

    Vector& Vector::operator=(Vector&& v)
    {
        LDELETE_ARRAY(x_);
        size_ = v.size_;
        x_ = v.x_;
        v.size_ = 0;
        v.x_ = NULL;
        return *this;
    }

    element_type dot(const Vector& v0, const Vector& v1)
    {
        LASSERT(v0.size() == v1.size());
        element_type d=0.0;
        for(s32 i=0; i<v0.size(); ++i){
            d += v0[i]*v1[i];
        }
        return d;
    }

    element_type dot(s32 size, const Vector& v0, const Vector& v1)
    {
        LASSERT(size <= v0.size());
        LASSERT(size <= v1.size());
        element_type d=0.0;
        for(s32 i=0; i<size; ++i){
            d += v0[i]*v1[i];
        }
        return d;
    }

    void print(const Vector& v)
    {
        for(s32 i=0; i<v.size(); ++i){
            printf("%f, ", v[i]);
        }
        printf("\n");
    }

    //--------------------------------------------------
    //---
    //--- VectorView
    //---
    //--------------------------------------------------
    VectorView::VectorView()
        :size_(0)
        ,x_(NULL)
    {
    }

    VectorView::VectorView(s32 n, element_type* x)
        :size_(n)
        ,x_(x)
    {
        LASSERT(0<=size_);
        LASSERT(NULL != x_);
    }

    VectorView::VectorView(VectorView&& v)
        :size_(v.size_)
        ,x_(v.x_)
    {
        v.size_ = 0;
        v.x_ = NULL;
    }

    VectorView::VectorView(Vector& v)
        :size_(v.size_)
        ,x_(v.x_)
    {
    }

    VectorView::~VectorView()
    {
        size_ = 0;
        x_ = NULL;
    }

    VectorView& VectorView::operator=(VectorView&& v)
    {
        size_ = v.size_;
        x_ = v.x_;
        v.size_ = 0;
        v.x_ = NULL;
        return *this;
    }


    //VectorView& VectorView::operator=(Vector& v)
    //{
    //    size_ = v.size_;
    //    x_ = v.x_;
    //    return *this;
    //}

    VectorView& VectorView::copy(const VectorView& v)
    {
        LASSERT(size_ == v.size());
        for(s32 i=0; i<size_; ++i){
            x_[i] = v[i];
        }
        return *this;
    }

    VectorView& VectorView::copy(const VectorStepView& v)
    {
        LASSERT(size_ == v.size());
        for(s32 i=0; i<size_; ++i){
            x_[i] = v[i];
        }
        return *this;
    }

    element_type VectorView::minimum() const
    {
        LASSERT(0<size_);
        element_type x = x_[0];
        for(s32 i=1; i<size_; ++i){
            if(x_[i]<x){
                x = x_[i];
            }
        }
        return x;
    }

    element_type VectorView::maximum() const
    {
        LASSERT(0<size_);
        element_type x = x_[0];
        for(s32 i=1; i<size_; ++i){
            if(x<x_[i]){
                x = x_[i];
            }
        }
        return x;
    }

    element_type innerproduct(const VectorView& v0, const VectorView& v1)
    {
        LASSERT(v0.size() == v1.size());
        element_type d=0.0;
        for(s32 i=0; i<v0.size(); ++i){
            d += v0[i]*v1[i];
        }
        return d;
    }

    element_type innerproduct(s32 offset, const VectorView& v0, const VectorView& v1)
    {
        LASSERT(v0.size() == v1.size());
        element_type norm=0.0;
        for(s32 i=offset; i<v0.size(); ++i){
            norm += v0[i]*v1[i];
        }
        return norm;
    }

    void print(const VectorView& v)
    {
        for(s32 i=0; i<v.size(); ++i){
            printf("%f, ", v[i]);
        }
        printf("\n");
    }

    //--------------------------------------------------
    //---
    //--- VectorStepView
    //---
    //--------------------------------------------------
    VectorStepView::VectorStepView()
        :size_(0)
        ,step_(0)
        ,x_(NULL)
    {
    }

    VectorStepView::VectorStepView(s32 n, s32 step, element_type* x)
        :size_(n)
        ,step_(step)
        ,x_(x)
    {
        LASSERT(0<=size_);
        LASSERT(0<=step);
        LASSERT(NULL != x_);
    }

    VectorStepView::VectorStepView(VectorStepView&& v)
        :size_(v.size_)
        ,step_(v.step_)
        ,x_(v.x_)
    {
        v.size_ = 0;
        v.step_ = 0;
        v.x_ = NULL;
    }

    VectorStepView::~VectorStepView()
    {
        size_ = 0;
        step_ = 0;
        x_ = NULL;
    }

    VectorStepView& VectorStepView::operator=(VectorStepView&& v)
    {
        size_ = v.size_;
        step_ = v.step_;
        x_ = v.x_;
        v.size_ = 0;
        v.step_ = 0;
        v.x_ = NULL;
        return *this;
    }

    VectorStepView& VectorStepView::copy(const VectorStepView& v)
    {
        LASSERT(size_ == v.size());
        for(s32 i=0; i<size_; ++i){
            (*this)[i] = v[i];
        }
        return *this;
    }

    VectorStepView& VectorStepView::copy(const VectorView& v)
    {
        LASSERT(size_ == v.size());
        for(s32 i=0; i<size_; ++i){
            (*this)[i] = v[i];
        }
        return *this;
    }

    element_type innerproduct(s32 offset, const VectorStepView& v0, const VectorStepView& v1)
    {
        LASSERT(v0.size() == v1.size());
        element_type norm=0.0;
        for(s32 i=offset; i<v0.size(); ++i){
            norm += v0[i]*v1[i];
        }
        return norm;
    }
}
