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
    element_type innerproduct(s32 size, const element_type* v0, const element_type* v1);

    class VectorStepView;

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
        inline element_type operator()(s32 index) const;
        inline element_type& operator()(s32 index);

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

    inline Vector::element_type Vector::operator()(s32 index) const
    {
        LASSERT(0<=index && index<size_);
        return x_[index];
    }

    inline Vector::element_type& Vector::operator()(s32 index)
    {
        LASSERT(0<=index && index<size_);
        return x_[index];
    }

    element_type dot(const Vector& v0, const Vector& v1);
    element_type dot(s32 size, const Vector& v0, const Vector& v1);
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
        inline element_type operator()(s32 index) const;
        inline element_type& operator()(s32 index);

        VectorView& operator=(const VectorView& v)
        {
            return copy(v);
        }
        VectorView& operator=(VectorView&& v);

        VectorView& operator=(const VectorStepView& v)
        {
            copy(v);
        }

        VectorView& copy(const VectorView& v);
        VectorView& copy(const VectorStepView& v);

        element_type minimum() const;
        element_type maximum() const;

    private:
        VectorView(const VectorView&) = delete;

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

    inline VectorView::element_type VectorView::operator()(s32 index) const
    {
        LASSERT(0<=index && index<size_);
        return x_[index];
    }

    inline VectorView::element_type& VectorView::operator()(s32 index)
    {
        LASSERT(0<=index && index<size_);
        return x_[index];
    }

    element_type innerproduct(const VectorView& v0, const VectorView& v1);
    element_type innerproduct(s32 offset, const VectorView& v0, const VectorView& v1);
    void print(const VectorView& v);


    //--------------------------------------------------
    //---
    //--- VectorStepView
    //---
    //--------------------------------------------------
    class VectorStepView
    {
    public:
        typedef lmath::element_type element_type;

        VectorStepView();
        VectorStepView(s32 n, s32 step, element_type* x);
        VectorStepView(VectorStepView&& v);
        ~VectorStepView();

        inline s32 size() const;
        inline element_type operator[](s32 index) const;
        inline element_type& operator[](s32 index);
        inline element_type operator()(s32 index) const;
        inline element_type& operator()(s32 index);

        VectorStepView& operator=(const VectorStepView& v)
        {
            return copy(v);
        }
        VectorStepView& operator=(VectorStepView&& v);

        VectorStepView& operator=(const VectorView& v)
        {
            return copy(v);
        }
        VectorStepView& copy(const VectorStepView& v);
        VectorStepView& copy(const VectorView& v);
    private:
        VectorStepView(const VectorStepView&) = delete;

        s32 size_;
        s32 step_;
        element_type* x_;
    };

    inline s32 VectorStepView::size() const
    {
        return size_;
    }

    inline VectorStepView::element_type VectorStepView::operator[](s32 index) const
    {
        LASSERT(0<=index && index<size_);
        return x_[step_*index];
    }

    inline VectorStepView::element_type& VectorStepView::operator[](s32 index)
    {
        LASSERT(0<=index && index<size_);
        return x_[step_*index];
    }

    inline VectorStepView::element_type VectorStepView::operator()(s32 index) const
    {
        LASSERT(0<=index && index<size_);
        return x_[step_*index];
    }

    inline VectorStepView::element_type& VectorStepView::operator()(s32 index)
    {
        LASSERT(0<=index && index<size_);
        return x_[step_*index];
    }

    element_type innerproduct(s32 offset, const VectorStepView& v0, const VectorStepView& v1);
}
#endif //INC_LMATH_VECTOR_H__