/**
@file Tridiagonalize.cs
@author t-sakai
@date 2017/06/04 create
*/

namespace lmath
{
    public static class Tridiagonalize
    {
        public static void tridiagonalize(int n, MatrixView m, VectorView d, int eoffset, VectorView e)
        {
            int n2 = n-2;
            int n1 = n-1;
            VectorView v,w;
            for(int k = 0; k<n2; ++k) {
                v = m.get(0, k);
                d[k] = v[k];
                e[k, eoffset] = Householder.householder(k+1, n1-k, v);
                if(-lmath.Math.Epsilon<e[k, eoffset] && e[k, eoffset]<lmath.Math.Epsilon) {
                    continue;
                }
                for(int i = k+1; i<n; ++i) {
                    double sum = 0.0;
                    for(int j = k+1; j<i; ++j) {
                        sum += m[i, j]*v[j];
                    }
                    for(int j = i; j<n; ++j) {
                        sum += m[j, i]*v[j];
                    }
                    d[i] = sum;
                }
                double t = Householder.innerproduct(k+1, n1-k, v, d)*0.5;
                for(int i = n1; k<i; --i) {
                    double p = v[i];
                    double q = d[i]-t*p;
                    d[i] = q;
                    for(int j = i; j<n; ++j) {
                        m[j, i] -= p*d[j] + q*v[j];
                    }
                }
            }//for(s32 k=0

            if(2<=n) {
                d[n2] = m[n2, n2];
                e[n2, eoffset] = m[n1, n2];
            }
            if(1<=n) {
                d[n1] = m[n1, n1];
            }
            for(int k = n1; 0<=k; --k) {
                v = m.get(0, k);
                if(k<n2) {
                    for(int i = k+1; i<n; ++i) {
                        w = m.get(0, i);
                        double t = Householder.innerproduct(k+1, n1-k, v, w);
                        for(int j = k+1; j<n; ++j) {
                            w[j] -= t*v[j];
                        }
                    }
                }
                for(int i = 0; i<n; ++i) {
                    v[i] = 0.0;
                }
                v[k] = 1.0;
            }
        }
    }
}
