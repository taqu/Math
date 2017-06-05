/**
@file Vector.cs
@author t-sakai
@date 2017/06/04 create
*/

namespace lmath
{
    //--------------------------------------------------
    //---
    //--- VectorView
    //---
    //--------------------------------------------------
    public struct VectorView
    {
        private int size_;
        private int offset_;
        private double[] x_;

        public VectorView(int size, int offset, double[] x)
        {
            size_ = size;
            offset_ = offset;
            x_ = x;
        }

        public int size()
        {
            return size_;
        }

        public double this[int index]
        {
            get { return x_[offset_+index]; }
            set { x_[offset_+index] = value; }
        }

        public double this[int index, int offset]
        {
            get { return x_[offset_+index+offset]; }
            set { x_[offset_+index+offset] = value; }
        }

        public void setSize(int size)
        {
            size_ = size;
        }

        public static double dot(VectorView v0, VectorView v1)
        {
            double d = 0.0;
            for(int i = 0; i<v0.size(); ++i) {
                d += v0[i]*v1[i];
            }
            return d;

        }

        public static double dot(VectorView v0, VectorView v1, int size)
        {
            double d = 0.0;
            for(int i = 0; i<size; ++i) {
                d += v0[i]*v1[i];
            }
            return d;
        }

        public static void print(VectorView v)
        {
            for(int i=0; i<v.size(); ++i){
                System.Console.Write(v[i]);
                System.Console.Write(", ");
            }
            System.Console.WriteLine();
        }
    };
}
