#include "Vector.h"
#include "Matrix.h"
#include "Householder.h"
#include "QR.h"
#include "SVD.h"
#include <stdio.h>
#include <random>
#include <iostream>

//#define USE_EIGEN
#ifdef USE_EIGEN
#include "../../eigen/Eigen/Householder"
#include "../../eigen/Eigen/QR"
#include "../../eigen/Eigen/SVD"


void copy(lmath::MatrixView& dst, Eigen::MatrixXd& src)
{
    for(lmath::s32 i=0; i<dst.rows(); ++i){
        for(lmath::s32 j=0; j<dst.cols(); ++j){
            dst(i,j) = src(i,j);
        }
    }
}
void transpose(lmath::MatrixView& dst, Eigen::MatrixXd& src)
{
    for(lmath::s32 i=0; i<dst.rows(); ++i){
        for(lmath::s32 j=0; j<dst.cols(); ++j){
            dst(j,i) = src(i,j);
        }
    }
}
#endif

void getRandom(lmath::MatrixView& m, std::mt19937& r)
{
    std::uniform_real_distribution<double> dist(0,10.0);
    for(int i=0; i<m.rows(); ++i){
        for(int j=0; j<m.cols(); ++j){
            m(i,j) = dist(r);
        }
    }
}

void getMSE(lmath::f64& mse, lmath::f64& minError, lmath::f64& maxError, lmath::MatrixView& m0, lmath::MatrixView& m1)
{
    minError = std::numeric_limits<lmath::f64>::max();
    maxError = 0.0;

    lmath::f64 t=0.0;
    for(int i=0; i<m0.rows(); ++i){
        for(int j=0; j<m0.cols(); ++j){
            lmath::f64 s = m0(i,j)-m1(i,j);
            t += s*s;
            s = lmath::absolute(s);
            if(s<minError){
                minError = s;
            }
            if(maxError<s){
                maxError = s;
            }
        }
    }
    mse = t/(m0.rows()*m0.cols());
}

bool equal(lmath::MatrixView& m0, lmath::MatrixView& m1, lmath::f64 epsilon)
{
    for(int i=0; i<m0.rows(); ++i){
        for(int j=0; j<m0.cols(); ++j){
            if(epsilon<lmath::absolute(m0(i,j)-m1(i,j))){
                return false;
            }
        }
    }
    return true;
}

bool checkBidiagonality(lmath::MatrixView& m, lmath::f64 epsilon)
{
    for(int j=0; j<m.cols()-1; ++j){
        for(int i=j+1; i<m.rows(); ++i){
            lmath::f64 t = lmath::absolute(m(i,j));
            if(epsilon<t){
                return false;
            }
        }
    }

    for(int j=2; j<m.cols(); ++j){
        int r = j-1;
        if(m.rows()<r){
            r = m.rows();
        }
        for(int i=0; i<r; ++i){
            lmath::f64 t = lmath::absolute(m(i,j));
            if(epsilon<t){
                return false;
            }
        }
    }
    return true;
}

int main(int argc, char** argv)
{
    std::random_device dev;
    std::mt19937 random(dev());

    {//Householder Transformation
        printf("Householder Transformation\n");
        printf("----------------------------------------\n");
        lmath::f64 b[10] = {10,9,8,7,6,5,4,3,2,1};
        lmath::f64 mbuffer0[10*10];
        lmath::f64 mbuffer1[10*10];
        lmath::f64 vbuffer0[10];
        lmath::f64 vbuffer1[10];

        lmath::VectorView b0(10, b);
        lmath::MatrixView m0(10,10,mbuffer0);
        lmath::MatrixView m1(10,10,mbuffer1);
        lmath::VectorView v0(10, vbuffer0);
        lmath::VectorView v1(10, vbuffer1);

        v0 = b0;
        lmath::householder_matrix(m0, &b0[0], 5, b0.size());
        lmath::print(m0);
        lmath::mul(m1, m0, m0);
        lmath::print(m1);
        lmath::mul(v1, m0, v0);
        lmath::print(v1);
    }

    {//Householder Transformation
        printf("----------------------------------------\n");
        lmath::f64 mbuffer0[4*4] = {4,1,-2,2, 1,2,0,1, -2,0,3,-2, 2,1,-2,-1};
        lmath::f64 mbuffer1[4*4];
        lmath::f64 mbuffer2[4*4];
        lmath::f64 vbuffer0[4];
        lmath::MatrixView m0(4,4,mbuffer0);
        lmath::MatrixView m1(4,4,mbuffer1);
        lmath::MatrixView m2(4,4,mbuffer2);
        lmath::VectorView r(4, vbuffer0);
        r.copy(m0.getRow(0));
        lmath::householder_matrix(m1, r, 1);
        m0 = lmath::mul(m2, m1, m0);
        m0 = lmath::mul(m2, m0, m1);
        printf("H\n");
        lmath::print(m1);
        printf("M\n");
        lmath::print(m0);
    }

    {//Bidiagonalization
        printf("Bidiagonalization\n");
        printf("----------------------------------------\n");
        lmath::f64 mbuffer0[4*4];
        lmath::f64 mbuffer1[4*4];
        lmath::f64 mbuffer2[4*4];
        lmath::f64 mbuffer3[4*4] = {1,3,2,1, 5,6,4,1, 7,8,9,1, 10,11,12,1};
        lmath::f64 mbuffer4[4*4];

        lmath::MatrixView U(4, 4, mbuffer0);
        lmath::MatrixView B(4, 4, mbuffer1);
        lmath::MatrixView V(4, 4, mbuffer2);
        lmath::MatrixView A(4, 4, mbuffer3);
        lmath::MatrixView T(4, 4, mbuffer4);

        lmath::bidiagonalization(U,B,V,A);

        std::cout << "U:\n";
        lmath::print(U);
        std::cout << "B:\n";
        lmath::print(B);
        std::cout << "V^T:\n";
        lmath::print(V);
        std::cout << "A:\n";
        lmath::print(A);
        A = lmath::mul(T, U, B);
        A = lmath::mul(T, A, V);
        std::cout << "U*B*V:\n";
        lmath::print(A);
    }
    {
        printf("----------------------------------------\n");
        lmath::f64 mbuffer0[4*4];
        lmath::f64 mbuffer1[4*3];
        lmath::f64 mbuffer2[3*3];
        lmath::f64 mbuffer3[4*3] = {1,3,2,5,6,4,7,8,9,10,11,12};
        lmath::f64 mbuffer4[4*3];

        lmath::MatrixView U(4, 4, mbuffer0);
        lmath::MatrixView B(4, 3, mbuffer1);
        lmath::MatrixView V(3, 3, mbuffer2);
        lmath::MatrixView A(4, 3, mbuffer3);
        lmath::MatrixView T(4, 3, mbuffer4);

        lmath::bidiagonalization(U,B,V,A);

        std::cout << "U:\n";
        lmath::print(U);
        std::cout << "B:\n";
        lmath::print(B);
        std::cout << "V^T:\n";
        lmath::print(V);
        std::cout << "A:\n";
        lmath::print(A);
        A = lmath::mul(T, U, B);
        A = lmath::mul(T, A, V);
        std::cout << "U*B*V:\n";
        lmath::print(A);
    }
    {
        printf("----------------------------------------\n");
        lmath::f64 mbuffer0[3*3];
        lmath::f64 mbuffer1[3*4];
        lmath::f64 mbuffer2[4*4];
        lmath::f64 mbuffer3[3*4] = {1,3,2,5,6,4,7,8,9,10,11,12};
        lmath::f64 mbuffer4[3*4];

        lmath::MatrixView U(3, 3, mbuffer0);
        lmath::MatrixView B(3, 4, mbuffer1);
        lmath::MatrixView V(4, 4, mbuffer2);
        lmath::MatrixView A(3, 4, mbuffer3);
        lmath::MatrixView T(3, 4, mbuffer4);

        lmath::bidiagonalization(U,B,V,A);

        std::cout << "U:\n";
        lmath::print(U);
        std::cout << "B:\n";
        lmath::print(B);
        std::cout << "V^T:\n";
        lmath::print(V);
        std::cout << "A:\n";
        lmath::print(A);
        A = lmath::mul(T, U, B);
        A = lmath::mul(T, A, V);
        std::cout << "U*B*V:\n";
        lmath::print(A);
    }
    {
        printf("----------------------------------------\n");
        const int Rows = 50;
        const int Cols = 50;
        lmath::f64 mbuffer0[Rows*Rows];
        lmath::f64 mbuffer1[Rows*Cols];
        lmath::f64 mbuffer2[Cols*Cols];
        lmath::f64 mbuffer3[Rows*Cols];
        lmath::f64 mbuffer4[Rows*Cols];

        lmath::MatrixView U(Rows, Rows, mbuffer0);
        lmath::MatrixView B(Rows, Cols, mbuffer1);
        lmath::MatrixView V(Cols, Cols, mbuffer2);
        lmath::MatrixView A(Rows, Cols, mbuffer3);
        lmath::MatrixView T(Rows, Cols, mbuffer4);

        for(int i = 0; i<2; ++i){
            getRandom(A, random);
            lmath::bidiagonalization(U, B, V, A);
            bool c = checkBidiagonality(B, 1.0e-13);
            B = lmath::mul(T, U, B);
            B = lmath::mul(T, B, V);
            std::cout << "[" << i << "] bidiagonal B " << c << "\n";

            if(!equal(A,B,std::numeric_limits<lmath::f64>::epsilon())){
                std::cout << "Error:\n";
                lmath::f64 mse, minError, maxError;
                getMSE(mse, minError, maxError, A,B);
                std::cout << "  mean squared error: " << mse << ", min: " << minError << ", max: " << maxError << std::endl;
            }
        }
    }

    {//Givens Rotation
        printf("Givens Rotation\n");
        printf("----------------------------------------\n");
        lmath::f64 f = 1.0;
        lmath::f64 g = 0.0;
        lmath::f64 c,s,r;
        lmath::rotGivens(c,s,r,f,g);
        printf("[f,g] = [%g, %g]\n [c,s,r] = [%g, %g, %g]\n", f, g, c, s, r);
        f = 1.0;
        g = 2.0;
        lmath::rotGivens(c,s,r,f,g);
        printf("[f,g] = [%g, %g]\n [c,s,r] = [%g, %g, %g]\n", f, g, c, s, r);
        lmath::f64 r0 = c*f+s*g;
        lmath::f64 r1 = -s*f+c*g;
        printf("[c s;-s c]*[f;g] = [r;0] : [%g;%g] = [%g;0]\n", r0, r1, r);
    }

    {//msweep,vsweep
        printf("msweep,vsweep\n");
        printf("----------------------------------------\n");
        lmath::f64 bbuffer[10*10] =
        {
            1, 11,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  2, 12,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  3, 13,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  4, 14,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  5, 15,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  6, 16,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  7, 17,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  8, 18,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  9, 19,
            0,  0,  0,  0,  0,  0,  0,  0,  0, 10,
        };
        lmath::f64 dbuffer[10] = {1,2,3,4,5,6,7,8,9,10};
        lmath::f64 ebuffer[9] = {11,12,13,14,15,16,17,18,19};
        lmath::MatrixView B(10, 10, bbuffer);
        lmath::VectorView d(10, dbuffer);
        lmath::VectorView e(9, ebuffer);
        lmath::msweep(B);
        printf("B:\n");
        lmath::print(B);

        lmath::vsweep(0, e.size(), d,e);
        printf("d:\n");
        lmath::print(d);
        printf("e:\n");
        lmath::print(e);
    }

    {//msweep2,vsweep2
        printf("msweep2,vsweep2\n");
        printf("----------------------------------------\n");
        lmath::f64 bbuffer[10*10] =
        {
            1, 11,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  2, 12,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  3, 13,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  4, 14,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  5, 15,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  6, 16,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  7, 17,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  8, 18,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  9, 19,
            0,  0,  0,  0,  0,  0,  0,  0,  0, 10,
        };
        lmath::f64 ubuffer[10*10];
        lmath::f64 vbuffer[10*10];
        lmath::f64 tbuffer[10*10];
        lmath::f64 dbuffer[10] = {1,2,3,4,5,6,7,8,9,10};
        lmath::f64 ebuffer[9] = {11,12,13,14,15,16,17,18,19};
        lmath::MatrixView B(10, 10, bbuffer);
        lmath::MatrixView U(10, 10, ubuffer);
        lmath::MatrixView V(10, 10, vbuffer);
        lmath::MatrixView T(10, 10, tbuffer);
        lmath::VectorView d(10, dbuffer);
        lmath::VectorView e(9, ebuffer);
        U.identity();
        V.identity();
        lmath::msweep2(1,9,U,B,V);
        printf("B:\n");
        lmath::print(B);
        printf("U:\n");
        lmath::print(U);
        printf("V:\n");
        lmath::print(V);
        printf("U*B*V:\n");
        B=lmath::mul(T,U,B);
        B=lmath::mul(T,B,V);
        lmath::print(B);

        lmath::vsweep(1,9,d,e);
        printf("d:\n");
        lmath::print(d);
        printf("e:\n");
        lmath::print(e);
    }

    {//Singular Value Decomposition
        printf("Singular Value Decomposition\n");
        printf("----------------------------------------\n");

        lmath::f64 bbuffer[10*10] =
        {
            1, 11,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  2, 12,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  3, 13,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  4, 14,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  5, 15,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  6, 16,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  7, 17,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  8, 18,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  9, 19,
            0,  0,  0,  0,  0,  0,  0,  0,  0, 10,
        };
        lmath::f64 abuffer[10*10];
        lmath::f64 apbuffer[10*10];
        lmath::f64 ubuffer[10*10];
        lmath::f64 utbuffer[10*10];
        lmath::f64 vbuffer[10*10];
        lmath::f64 vtbuffer[10*10];
        lmath::f64 tbuffer[10*10];
        lmath::f64 dbuffer[10] = {1,2,3,4,5,6,7,8,9,10};
        lmath::f64 ebuffer[9] = {11,12,13,14,15,16,17,18,19};
        lmath::MatrixView B(10, 10, bbuffer);
        lmath::MatrixView A(10, 10, abuffer);
        lmath::MatrixView Aplus(10, 10, apbuffer);
        lmath::MatrixView U(10, 10, ubuffer);
        lmath::MatrixView Ut(10, 10, utbuffer);
        lmath::MatrixView V(10, 10, vbuffer);
        lmath::MatrixView Vt(10, 10, vtbuffer);
        lmath::MatrixView T(10, 10, tbuffer);
        lmath::VectorView d(10, dbuffer);
        lmath::VectorView e(9, ebuffer);

        lmath::s32 iterations=0;
        iterations = lmath::svd(d,e);
        printf("svd result: %d\n", iterations);
        printf("d:\n");
        lmath::print(d);
        printf("e:\n");
        lmath::print(e);

        iterations = lmath::svd(U,B,V);
        printf("svd result: %d\n", iterations);
        A = B;
        U.transpose(Ut);
        V.transpose(Vt);
        printf("B:\n");
        lmath::print(B);
        printf("U:\n");
        lmath::print(U);
        printf("V:\n");
        lmath::print(V);
        printf("U*B*V:\n");
        B=lmath::mul(T,U,B);
        B=lmath::mul(T,B,V);
        lmath::print(B);

        lmath::pseudoInverse(Aplus, U, A, V);
        printf("B+:\n");
        lmath::print(Aplus);

        printf("U*Ut:\n");
        Ut = lmath::mul(T, U, Ut);
        lmath::print(Ut);
        printf("V*Vt:\n");
        Vt = lmath::mul(T, V, Vt);
        lmath::print(Vt);
    }

     {//Singular Value Decomposition
        printf("Solve Linear Equation\n");
        printf("----------------------------------------\n");

        lmath::f64 bbuffer[3*2] =
        {
            1,2,
            4,5,
            7,8,
        };
        lmath::f64 vbbuffer[3] =
        {
            3,
            6,
            9,
        };

        lmath::f64 bpbuffer[3*2];
        lmath::f64 ubuffer[3*3];
        lmath::f64 vbuffer[2*2];
        lmath::f64 tbuffer[3*2];
        lmath::f64 tvbuffer[2];
        lmath::MatrixView B(3, 2, bbuffer);
        lmath::MatrixView Bplus(2, 3, bpbuffer);
        lmath::MatrixView U(3, 3, ubuffer);
        lmath::MatrixView V(2, 2, vbuffer);
        lmath::MatrixView T(3, 2, tbuffer);
        lmath::VectorView b(3,vbbuffer);
        lmath::VectorView t(2, tvbuffer);


        lmath::s32 iterations=0;
        iterations = lmath::svd(U,B,V);
        printf("svd result: %d\n", iterations);
        lmath::pseudoInverse(Bplus, U, B, V);
        printf("B+:\n");
        lmath::print(Bplus);
        lmath::mul(t, Bplus, b);
        printf("x:\n");
        lmath::print(t);
    }
    return 0;
}
