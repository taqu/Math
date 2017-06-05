/**
@file Matrix.cs
@author t-sakai
@date 2017/06/04 create
*/

namespace lmath
{
    //--------------------------------------------------
    //---
    //--- MatrixView
    //---
    //--------------------------------------------------
    public struct MatrixView
    {
        private int cols_;
        private int rows_;
        private int offset_;
        private double[] x_;

        public MatrixView(int cols, int rows, int offset, double[] x)
        {
            cols_ = cols;
            rows_ = rows;
            offset_ = offset;
            x_ = x;
        }

        public int cols()
        {
            return cols_;
        }
        public int rows()
        {
            return rows_;
        }
        public double this[int c, int r]
        {
            get{ return x_[offset_ + r*cols_ + c]; }
            set{ x_[offset_ + r*cols_ + c] = value; }
        }

        public VectorView get(int c, int r)
        {
            int index = r*cols_ + c;
            int size = cols_*rows_-index;
            return new VectorView(size, index+offset_, x_);
        }

        public void setSize(int cols, int rows)
        {
            cols_ = cols;
            rows_ = rows;
        }

        public void setIdentity()
        {
            System.Array.Clear(x_, offset_, cols_ * rows_);
            int n = (cols_<rows_)? cols_ : rows_;
            for(int i=0; i<n; ++i){
                this[i,i] = 1.0;
            }
        }

        public void transpose()
        {
            for(int i = 0; i < rows_; ++i){
                for(int j = 0; j < cols_; ++j){
                    int i0 = offset_ + j*cols_+i;
                    int i1 = offset_ + i*cols_+j;
                    lmath.Math.swap(ref x_[i0], ref x_[i1]);
                }
            }
        }

        public void transpose(MatrixView dst)
        {
            for(int i = 0; i < rows_; ++i){
                for(int j = 0; j < cols_; ++j){
                    dst[i,j] = x_[offset_ + i*cols_ + j];
                }
            }
        }

        public static void mul(MatrixView dst, MatrixView m0, MatrixView m1)
        {
            for(int i=0; i<m0.rows(); ++i){
                for(int j=0; j<m1.cols(); ++j){
                    double t = 0.0;
                    for(int k=0; k < m0.cols(); ++k){
                        t += m0[k, i] * m1[j, k];
                    }
                    dst[j, i] = t;
                }
            }
        }

        public static void print(MatrixView m)
        {
            for(int i=0 ; i < m.rows(); ++i){
                for(int j=0; j < m.cols(); ++j){
                    System.Console.Write(m[j, i]);
                    System.Console.Write(", ");
                }
                System.Console.WriteLine();
            }
        }
    };
}
