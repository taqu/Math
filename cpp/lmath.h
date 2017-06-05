#ifndef INC_LMATH_LMATH_H__
#define INC_LMATH_LMATH_H__
/**
@file lmath.h
@author t-sakai
@date 2017/06/02 create
*/
#include <cassert>
#include <cmath>
#include <type_traits>

namespace lmath
{
    typedef __int32 s32;
    typedef __int64 s64;

    typedef unsigned __int32 u32;
    typedef unsigned __int64 u64;

    typedef float f32;
    typedef double f64;

#ifndef NULL
#define NULL (0)
#endif
#define LASSERT(exp) assert(exp)
#define LNEW new
#define LDELETE(ptr) delete (ptr); (ptr)=NULL
#define LDELETE_ARRAY(ptr) delete[] (ptr); (ptr)=NULL

#define LMATH_EPSILON (1.0e-6f)
#define LMATH_TRUNCATE_EPSILON (1.0e-6f)

    union UnionU32F32
    {
        u32 u32_;
        f32 f32_;
    };

    union UnionU64F64
    {
        u64 u64_;
        f64 f64_;
    };

    template<class T>
    inline T absolute(T x)
    {
        return (x<0)? -x : x;
    }

    template<>
    inline s32 absolute<s32>(s32 x)
    {
        x &= 0x7FFFFFFFU;
        return x;
    }

    template<>
    inline u32 absolute<u32>(u32 x)
    {
        return x;
    }

    template<>
    inline f32 absolute<f32>(f32 x)
    {
        UnionU32F32 u;
        u.f32_ = x;
        u.u32_ &= 0x7FFFFFFFU;
        return u.f32_;
    }

    template<>
    inline f64 absolute<f64>(f64 x)
    {
        UnionU64F64 u;
        u.f64_ = x;
        u.u64_ &= 0x7FFFFFFFFFFFFFFFU;
        return u.f64_;
    }

    typedef f32 element_type;

    inline f32 sqrt(f32 x)
    {
        return ::sqrtf(x);
    }

    inline f64 sqrt(f64 x)
    {
        return ::sqrt(x);
    }

    using std::move;

    template<class T>
    inline void swap(T& x0, T& x1)
    {
        T t(lmath::move(x0));
        x0 = lmath::move(x1);
        x1 = lmath::move(t);
    }
}
#endif //INC_LMATH_LMATH_H__
