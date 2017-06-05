/**
@file QR.cs
@author t-sakai
@date 2017/06/04 create
*/

namespace lmath
{
    public static class QR
    {
        public static int eigen(int n, MatrixView m, VectorView d, VectorView e, double epsilon=lmath.Math.Epsilon, int maxIteration=10)
        {
            Tridiagonalize.tridiagonalize(n, m, d, 1, e);
            e[0] = 0.0;
            for(int h = n-1; 0<h; --h) {
                int j = h;
                while(0<j && epsilon * (Math.absolute(d[j-1])+Math.absolute(d[j])) < Math.absolute(e[j])) {
                    --j;
                }
                if(j==h) {
                    continue;
                }
                int iteration = 0;
                double e0, e1;
                do {
                    if(maxIteration<++iteration) {
                        return -1;
                    }
                    double w = (d[h-1] - d[h])*0.5;
                    double t = e[h]*e[h];
                    double s = System.Math.Sqrt(w*w+t);
                    if(w<0.0) {
                        s = -s;
                    }
                    double x = d[j] - d[h] + t/(w+s);
                    double y = e[j+1];
                    for(int k = j; k<h; ++k) {
                        double c;
                        if(Math.absolute(y)<=Math.absolute(x)) {
                            t = -y/x;
                            c = 1.0/System.Math.Sqrt(t*t+1.0);
                            s = t*c;
                        } else {
                            t = -x/y;
                            s = 1.0/System.Math.Sqrt(t*t+1.0);
                            c = t*s;
                        }
                        w = d[k] - d[k+1];
                        t = (w*s + 2.0*c*e[k+1])*s;
                        d[k] -= t;
                        d[k+1] += t;
                        if(j<k) {
                            e[k] = c*e[k] - s*y;
                        }
                        e[k+1] += s*(c*w - 2.0*s*e[k+1]);

                        //固有ベクトルを計算
                        for(int i = 0; i<n; ++i) {
                            x = m[i, k];
                            y = m[i, k+1];
                            m[i, k] = c*x - s*y;
                            m[i, k+1] = s*x + c*y;
                        }
                        if(k<(h-1)) {
                            x = e[k+1];
                            y = -s*e[k+2];
                            e[k+2] *= c;
                        }
                    }//for(int k=j;

                    e0 = epsilon*(Math.absolute(d[h-1])+Math.absolute(d[h]));
                    e1 = Math.absolute(e[h]);
                } while(e0<=e1);
            }

            //固有値, 固有ベクトルを降順に整列
            for(int k = 0; k<n-1; ++k) {
                int h = k;
                double t = d[h];
                for(int i = k+1; i<n; ++i) {
                    if(t<d[i]) {
                        h=i;
                        t = d[h];
                    }
                    d[h] = d[k];
                    d[k] = t;
                    for(int j = 0; j<m.cols(); ++j) {
                        t = m[j, h];
                        m[j, h] = m[j, k];
                        m[j, k] = t;
                    }
                }
            }
            return 0;
        }
    }
}
