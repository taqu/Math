/**
@file Math.cs
@author t-sakai
@date 2017/06/04 create
*/

namespace lmath
{
    public static class Math
    {
        public const double Epsilon = 1.0e-13;
        public const double TruncateEpsilon = 1.0e-10;

        public static void swap(ref double x0, ref double x1)
        {
            double tmp = x0;
            x0 = x1;
            x1 = tmp;
        }

        public static double absolute(double x)
        {
            return x<0.0? -x : x;
        }
    }
}
