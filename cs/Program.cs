using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Math
{
    class Program
    {
        static void Main(string[] args)
        {
            double[] buffer = new double[1024];
            {
                lmath.VectorView v = new lmath.VectorView(3, 0, buffer);
                v[0] = 1.0; v[1] = 1.0; v[2] = 1.0;
                double x0 = lmath.Householder.householder(0, v.size(), v);
                System.Console.Write(x0);
                System.Console.Write(" 0 0\n");
                System.Console.WriteLine();
            }
            {
                lmath.MatrixView original = new lmath.MatrixView(4, 4, 0, buffer);
                original[0, 0] = 3;
                original[1, 0] = 2;
                original[2, 0] = 4;
                original[3, 0] = 1;

                original[0, 1] = 2;
                original[1, 1] = 2;
                original[2, 1] = 4;
                original[3, 1] = 5;

                original[0, 2] = 4;
                original[1, 2] = 4;
                original[2, 2] = 1;
                original[3, 2] = 3;

                original[0, 3] = 1;
                original[1, 3] = 5;
                original[2, 3] = 3;
                original[3, 3] = 2;

                lmath.VectorView d = new lmath.VectorView(4, 16, buffer);
                lmath.VectorView e = new lmath.VectorView(4, 20, buffer);

                lmath.MatrixView.print(original);
                System.Console.WriteLine();
                lmath.Tridiagonalize.tridiagonalize(original.cols(), original, d, 0, e);
                lmath.VectorView.print(d);
                System.Console.WriteLine();
                lmath.VectorView.print(e);
                System.Console.WriteLine();
            }

            {
                lmath.MatrixView o = new lmath.MatrixView(3, 3, 0, buffer);
                lmath.MatrixView m = new lmath.MatrixView(3, 3, 9, buffer);
                lmath.VectorView d = new lmath.VectorView(3, 18, buffer);
                lmath.VectorView e = new lmath.VectorView(3, 21, buffer);
                o[0, 0] = m[0, 0] = 0.0f;
                o[1, 0] = m[1, 0] = 2.0f;
                o[2, 0] = m[2, 0] = 2.0f;

                o[0, 1] = m[0, 1] = 2.0f;
                o[1, 1] = m[1, 1] = 1.0f;
                o[2, 1] = m[2, 1] = 0.0f;

                o[0, 2] = m[0, 2] = 2.0f;
                o[1, 2] = m[1, 2] = 0.0f;
                o[2, 2] = m[2, 2] = -1.0f;

                lmath.MatrixView.print(m);
                System.Console.WriteLine();
                lmath.QR.eigen(m.cols(), m, d, e);
                lmath.MatrixView.print(m);
                System.Console.WriteLine();
                lmath.VectorView.print(d);
                System.Console.WriteLine();
                lmath.VectorView.print(e);
                System.Console.WriteLine();

                double x0, x1, x2;
                x0 = m[0, 0]*o[0, 0] + m[1, 0]*o[0, 1] + m[2, 0]*o[0, 2];
                x1 = m[0, 0]*o[1, 0] + m[1, 0]*o[1, 1] + m[2, 0]*o[1, 2];
                x2 = m[0, 0]*o[2, 0] + m[1, 0]*o[2, 1] + m[2, 0]*o[2, 2];
                System.Console.WriteLine(string.Format("{0}, {1}, {2}", x0, x1, x2));
                System.Console.WriteLine();

                x0 = m[0, 1]*o[0, 0] + m[1, 1]*o[0, 1] + m[2, 1]*o[0, 2];
                x1 = m[0, 1]*o[1, 0] + m[1, 1]*o[1, 1] + m[2, 1]*o[1, 2];
                x2 = m[0, 1]*o[2, 0] + m[1, 1]*o[2, 1] + m[2, 1]*o[2, 2];
                System.Console.WriteLine(string.Format("{0}, {1}, {2}", x0, x1, x2));
                System.Console.WriteLine();

                x0 = m[0, 2]*o[0, 0] + m[1, 2]*o[0, 1] + m[2, 2]*o[0, 2];
                x1 = m[0, 2]*o[1, 0] + m[1, 2]*o[1, 1] + m[2, 2]*o[1, 2];
                x2 = m[0, 2]*o[2, 0] + m[1, 2]*o[2, 1] + m[2, 2]*o[2, 2];
                System.Console.WriteLine(string.Format("{0}, {1}, {2}", x0, x1, x2));
                System.Console.WriteLine();
            }

            {
                const int NumSamples = 16;
                System.Random random = new System.Random();

                double a = random.NextDouble() * 10.0 + 1.0;
                double b = random.NextDouble() * 10.0 + 1.0;
                double[] x = new double[NumSamples];
                double[] y = new double[NumSamples];
                for(int i = 0; i<NumSamples; ++i) {
                    x[i] = random.NextDouble() * 100.0;
                    y[i] = a* x[i] + b + random.NextDouble()*0.2-0.1;
                }

                lmath.SVD svd = new lmath.SVD(2, NumSamples);
                lmath.MatrixView ma = svd.a();
                lmath.VectorView vb = svd.b();
                for(int i = 0; i<NumSamples; ++i) {
                    ma[0, i] = x[i];
                    ma[1, i] = 1.0;
                    vb[i] = y[i];
                }
                svd.solve();
                System.Console.WriteLine(string.Format("y={0} x + {1} \n", a, b));
                System.Console.Write("a,b = ");
                lmath.VectorView.print(svd.x());
                System.Console.WriteLine();
                lmath.VectorView vx = svd.x();
                double error = 0.0;
                for(int i=0; i<NumSamples; ++i){
                    double e = vx[0]*x[i] + vx[1] - y[i];
                    error += e*e;
                }
                System.Console.WriteLine(string.Format("error = {0}\n", error));
            }
        }
    }
}
