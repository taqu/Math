/**
@file SVD.cs
@author t-sakai
@date 2017/06/04 create
*/

namespace lmath
{
    //-----------------------------------------
    //---
    //--- SVD
    //---
    //-----------------------------------------
    public class SVD
    {
        private int cols_;
        private int rows_;
        private double[] buffer_;

        private MatrixView a_;
        private MatrixView at_;
        private MatrixView ap_;
        private MatrixView ut_;
        private MatrixView v_;
        private MatrixView sp_;
        private VectorView sigma_;
        private VectorView e_;

        private VectorView x_;
        private VectorView b_;

        public SVD(int cols, int rows)
        {
            reset(cols, rows);
        }

        public MatrixView a()
        {
            return a_;
        }

        public MatrixView ut()
        {
            return ut_;
        }


        public MatrixView v()
        {
            return v_;
        }

        public VectorView sigma()
        {
            return sigma_;
        }


        public VectorView x()
        {
            return x_;
        }

        public VectorView b()
        {
            return b_;
        }

        public void reset(int cols, int rows)
        {
            if(cols_ == cols && rows_ == rows) {
                return;
            }
            int size0 = cols*rows;
            int size1 = rows*rows;
            int size2 = cols*cols;

            int n = (rows<cols) ? cols : rows;
            if(size0 != (cols_*rows_)) {
                int total = size0*4 + size1 + size2 + n*2 + cols + rows;
                buffer_ = new double[total];
            }
            cols_ = cols;
            rows_ = rows;

            int offset = 0;

            a_ = new MatrixView(cols_, rows_, offset, buffer_);
            offset += size0;

            at_ = new MatrixView(rows_, cols_, offset, buffer_);
            offset += size0;

            ap_ = new MatrixView(rows_, rows_, offset, buffer_);
            offset += size0;

            ut_ = new MatrixView(rows_, rows_, offset, buffer_);
            offset += size1;

            v_ = new MatrixView(cols_, cols_, offset, buffer_);
            offset += size2;

            sp_ = new MatrixView(rows_, cols_, offset, buffer_);
            offset += size0;

            sigma_ = new VectorView(n, offset, buffer_);
            offset += n;

            e_ = new VectorView(n, offset, buffer_);
            offset += n;

            x_ = new VectorView(cols_, offset, buffer_);
            offset += cols_;

            b_ = new VectorView(rows_, offset, buffer_);
            offset += rows_;
        }

        public int solve(double epsilon = Math.Epsilon, int maxIteration = 10, double truncate = Math.TruncateEpsilon)
        {
            a_.transpose(at_);

            int ret = 0;
            MatrixView.mul(ut_, a_, at_);
            if(0 != QR.eigen(ut_.cols(), ut_, sigma_, e_, epsilon, maxIteration)) {
                ret |= 0x01;
            }

            MatrixView.mul(v_, at_, a_);
            if(0 != QR.eigen(v_.cols(), v_, sigma_, e_, epsilon, maxIteration)) {
                ret |= 0x02;
            }
            v_.transpose();

            sp_.setIdentity();
            for(int i = 0; i<cols_; ++i) {
                if(truncate<sigma_[i]) {
                    sp_[i, i] = 1.0f/System.Math.Sqrt(sigma_[i]);
                }else {
                    sp_[i, i] = 0.0;
                }
            }
            MatrixView.mul(at_, sp_, ut_);
            MatrixView.mul(ap_, v_, at_);
            for(int i = 0; i<cols_; ++i) {
                x_[i] = 0.0;
                for(int j = 0; j<rows_; ++j) {
                    x_[i] += ap_[j, i] * b_[j];
                }
            }
            return ret;
        }
    }
}
