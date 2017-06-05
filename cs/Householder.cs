/**
@file Householder.cs
@author t-sakai
@date 2017/06/04 create
*/

namespace lmath
{
    public static class Householder
    {
        public static double innerproduct(int offset, int size, VectorView v)
        {
            double d = 0.0;
            for(int i = 0; i<size; ++i) {
                double x = v[i, offset];
                d += x*x;
            }
            return d;
        }

        public static double innerproduct(int offset, int size, VectorView v0, VectorView v1)
        {
            double d = 0.0;
            for(int i = 0; i<size; ++i) {
                d += v0[i, offset]*v1[i, offset];
            }
            return d;
        }

        public static double householder(int offset, int size, VectorView v)
        {
            double norm = System.Math.Sqrt(innerproduct(offset, size, v));
            if(Math.Epsilon<=norm) {
                if(v[0+offset]<0.0f) {
                    norm = -norm;
                }
                v[0+offset] += norm;
                double weight = 1.0/System.Math.Sqrt(norm*v[0+offset]);
                for(int i = 0; i<size; ++i) {
                    v[i+offset] *= weight;
                }
            }
            return -norm;
        }
    }
}
