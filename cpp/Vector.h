#ifndef INC_LMATH_VECTOR_H__
#define INC_LMATH_VECTOR_H__
/**
@file Vector.h
@author t-sakai
@date 2017/06/02 create
*/
#include "lmath.h"

namespace lmath
{
    element_type innerproduct(s32 size, element_type* v0, element_type* v1);

    //--------------------------------------------------
    //---
    //--- Vector
    //---
    //--------------------------------------------------
    class Vector
    {
    public:
        typedef lmath::element_type element_type;

        Vector();
        explicit Vector(s32 n);
        Vector(const Vector& v);
        Vector(Vector&& v);
        ~Vector();

        inline s32 size() const;
        inline element_type operator[](s32 index) const;
        inline element_type& operator[](s32 index);

        Vector& operator=(const Vector& v);
        Vector& operator=(Vector&& v);

    private:
        friend class VectorView;

        s32 size_;
        element_type* x_;
    };

    inline s32 Vector::size() const
    {
        return size_;
    }

    inline Vector::element_type Vector::operator[](s32 index) const
    {
        LASSERT(0<=index && index<size_);
        return x_[index];
    }

    inline Vector::element_type& Vector::operator[](s32 index)
    {
        LASSERT(0<=index && index<size_);
        return x_[index];
    }

    element_type dot(const Vector& v0, const Vector& v1);
    element_type dot(const Vector& v0, const Vector& v1, s32 size);
    void print(const Vector& v);

    //--------------------------------------------------
    //---
    //--- VectorView
    //---
    //--------------------------------------------------
    class VectorView
    {
    public:
        typedef lmath::element_type element_type;

        VectorView();
        VectorView(s32 n, element_type* x);
        VectorView(VectorView&& v);
        VectorView(Vector& v);
        ~VectorView();

        inline s32 size() const;
        inline element_type operator[](s32 index) const;
        inline element_type& operator[](s32 index);

        VectorView& operator=(VectorView&& v);
        VectorView& operator=(Vector& v);

        void setSize(s32 size);
    private:
        VectorView(const VectorView&) = delete;
        VectorView& operator=(const VectorView&) = delete;

        s32 size_;
        element_type* x_;
    };

    inline s32 VectorView::size() const
    {
        return size_;
    }

    inline VectorView::element_type VectorView::operator[](s32 index) const
    {
        LASSERT(0<=index && index<size_);
        return x_[index];
    }

    inline VectorView::element_type& VectorView::operator[](s32 index)
    {
        LASSERT(0<=index && index<size_);
        return x_[index];
    }

    element_type dot(const VectorView& v0, const VectorView& v1);
    element_type dot(const VectorView& v0, const VectorView& v1, s32 size);
    void print(const VectorView& v);
}
#endif //INC_LMATH_VECTOR_H__