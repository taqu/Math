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
    element_type innerproduct(s32 size, element_type* v0, element_type* v1)
    {
        LASSERT(NULL != v0);
        LASSERT(NULL != v1);
        element_type d=0.0;
        for(s32 i=0; i<size; ++i){
            d += v0[i]*v1[i];
        }
        return d;
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
        LDELETE_ARRAY(x_);
        size_ = v.size_;
        x_ = LNEW element_type[size_];
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

    element_type dot(const Vector& v0, const Vector& v1, s32 size)
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
    }

    VectorView& VectorView::operator=(VectorView&& v)
    {
        LDELETE_ARRAY(x_);
        size_ = v.size_;
        x_ = v.x_;
        v.size_ = 0;
        v.x_ = NULL;
        return *this;
    }


    VectorView& VectorView::operator=(Vector& v)
    {
        size_ = v.size_;
        x_ = v.x_;
        return *this;
    }

    void VectorView::setSize(s32 size)
    {
        LASSERT(size<=size_);
        size_ = size;
    }

    element_type dot(const VectorView& v0, const VectorView& v1)
    {
        LASSERT(v0.size() == v1.size());
        element_type d=0.0;
        for(s32 i=0; i<v0.size(); ++i){
            d += v0[i]*v1[i];
        }
        return d;
    }

    element_type dot(const VectorView& v0, const VectorView& v1, s32 size)
    {
        LASSERT(size <= v0.size());
        LASSERT(size <= v1.size());
        element_type d=0.0;
        for(s32 i=0; i<size; ++i){
            d += v0[i]*v1[i];
        }
        return d;
    }

    void print(const VectorView& v)
    {
        for(s32 i=0; i<v.size(); ++i){
            printf("%f, ", v[i]);
        }
        printf("\n");
    }
}
