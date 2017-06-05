#include "Vector.h"
#include "Matrix.h"
#include "Householder.h"
#include "Tridiagonalize.h"
#include "QR.h"
#include "SVD.h"
#include <stdio.h>
#include <random>

int main(int argc, char** argv)
{
    lmath::Vector v(3);
    v[0] = 1.0f; v[1] = 1.0f; v[2] = 1.0f;

    lmath::f32 x0 = lmath::householder(v);

    printf("%f 0 0\n", x0);
    printf("\n");

    {
        lmath::Matrix original(4, 4);

        original(0, 0) = 3;
        original(1, 0) = 2;
        original(2, 0) = 4;
        original(3, 0) = 1;

        original(0, 1) = 2;
        original(1, 1) = 2;
        original(2, 1) = 4;
        original(3, 1) = 5;

        original(0, 2) = 4;
        original(1, 2) = 4;
        original(2, 2) = 1;
        original(3, 2) = 3;

        original(0, 3) = 1;
        original(1, 3) = 5;
        original(2, 3) = 3;
        original(3, 3) = 2;

        lmath::print(original);
        printf("\n");
        lmath::tridiagonalize(original);
        lmath::print(original);
        printf("\n");

    }

    {
        lmath::Matrix original(4, 4);
        lmath::Vector d(4);
        lmath::Vector e(4);

        original(0, 0) = 3;
        original(1, 0) = 2;
        original(2, 0) = 4;
        original(3, 0) = 1;

        original(0, 1) = 2;
        original(1, 1) = 2;
        original(2, 1) = 4;
        original(3, 1) = 5;

        original(0, 2) = 4;
        original(1, 2) = 4;
        original(2, 2) = 1;
        original(3, 2) = 3;

        original(0, 3) = 1;
        original(1, 3) = 5;
        original(2, 3) = 3;
        original(3, 3) = 2;

        lmath::print(original);
        printf("\n");
        lmath::MatrixView mview(original);
        lmath::VectorView dview(d);
        lmath::VectorView eview(e);

        lmath::tridiagonalize(4, mview, &dview[0], &eview[0]);
        lmath::print(dview);
        printf("\n");
        lmath::print(eview);
        printf("\n");
        lmath::print(mview);
        printf("\n");
    }

    {
        lmath::Matrix m(3, 3);
        m(0,0) = 0.0f;
        m(1,0) = 2.0f;
        m(2,0) = 2.0f;

        m(0,1) = 2.0f;
        m(1,1) = 1.0f;
        m(2,1) = 0.0f;

        m(0,2) = 2.0f;
        m(1,2) = 0.0f;
        m(2,2) = -1.0f;

        lmath::print(m);
        printf("\n");
        lmath::qr_algorithm(m);
        lmath::print(m);
        printf("\n");
    }

    {
        lmath::Matrix o(3, 3);
        lmath::Matrix m(3, 3);
        lmath::Vector d(3);
        lmath::Vector e(3);
        o(0,0) = m(0,0) = 0.0f;
        o(1,0) = m(1,0) = 2.0f;
        o(2,0) = m(2,0) = 2.0f;

        o(0,1) = m(0,1) = 2.0f;
        o(1,1) = m(1,1) = 1.0f;
        o(2,1) = m(2,1) = 0.0f;

        o(0,2) = m(0,2) = 2.0f;
        o(1,2) = m(1,2) = 0.0f;
        o(2,2) = m(2,2) = -1.0f;

        lmath::MatrixView mview(m);
        lmath::VectorView dview(d);
        lmath::VectorView eview(e);

        lmath::print(mview);
        printf("\n");
        lmath::eigen(mview.cols(), mview, dview, eview);
        lmath::print(mview);
        printf("\n");
        lmath::print(dview);
        printf("\n");
        lmath::print(eview);
        printf("\n");

        lmath::element_type x0,x1,x2;
        x0 = mview(0,0)*o(0,0) + mview(1,0)*o(0,1) + mview(2,0)*o(0,2);
        x1 = mview(0,0)*o(1,0) + mview(1,0)*o(1,1) + mview(2,0)*o(1,2);
        x2 = mview(0,0)*o(2,0) + mview(1,0)*o(2,1) + mview(2,0)*o(2,2);
        printf("%f, %f, %f\n\n", x0, x1, x2);

        x0 = mview(0,1)*o(0,0) + mview(1,1)*o(0,1) + mview(2,1)*o(0,2);
        x1 = mview(0,1)*o(1,0) + mview(1,1)*o(1,1) + mview(2,1)*o(1,2);
        x2 = mview(0,1)*o(2,0) + mview(1,1)*o(2,1) + mview(2,1)*o(2,2);
        printf("%f, %f, %f\n\n", x0, x1, x2);

        x0 = mview(0,2)*o(0,0) + mview(1,2)*o(0,1) + mview(2,2)*o(0,2);
        x1 = mview(0,2)*o(1,0) + mview(1,2)*o(1,1) + mview(2,2)*o(1,2);
        x2 = mview(0,2)*o(2,0) + mview(1,2)*o(2,1) + mview(2,2)*o(2,2);
        printf("%f, %f, %f\n\n", x0, x1, x2);
    }

    {
        lmath::Matrix tmp_m(3, 4);
        tmp_m(0, 0) = 0.411303; tmp_m(1, 0) = 0.409863; tmp_m(2, 0) = 0.278012;
        tmp_m(0, 1) = 0.246353; tmp_m(1, 1) = 0.098114; tmp_m(2, 1) = 0.907056;
        tmp_m(0, 2) = 0.569562; tmp_m(1, 2) = 0.813442; tmp_m(2, 2) = 0.163729;
        tmp_m(0, 3) = 0.757962; tmp_m(1, 3) = 0.873716; tmp_m(2, 3) = 0.260235;
        lmath::MatrixView m(tmp_m);

        lmath::Matrix tmp_u(4, 4);
        lmath::Matrix tmp_vt(3, 3);
        lmath::Vector tmp_sigma(4);

        lmath::MatrixView u(tmp_u);
        lmath::MatrixView vt(tmp_vt);
        lmath::VectorView sigma(tmp_sigma);
        lmath::svd(u, sigma, vt, m);

        printf("U\n");
        print(u);
        printf("\n");

        printf("Vt\n");
        print(vt);
        printf("\n");

        printf("sigma\n");
        print(sigma);
        printf("\n");
    }

    {
        static const lmath::s32 NumSamples = 16;
        std::random_device device;
        std::mt19937 random(device());
        std::uniform_real_distribution<lmath::f32> distX(0.0f, 100.0f);
        std::uniform_real_distribution<lmath::f32> distA(1.0f, 10.0f);
        std::uniform_real_distribution<lmath::f32> distError(-0.1f, 0.1f);

        lmath::f32 a = distA(random);
        lmath::f32 b = distA(random);
        lmath::f32 x[NumSamples];
        lmath::f32 y[NumSamples];
        for(lmath::s32 i=0; i<NumSamples; ++i){
            x[i] = distX(random);
            y[i] = a*x[i] + b + distError(random);
        }

        lmath::SVD svd(2, NumSamples);
        lmath::MatrixView& ma = svd.a();
        lmath::VectorView& vb = svd.b();
        for(lmath::s32 i=0; i<NumSamples; ++i){
            ma(0,i) = x[i];
            ma(1,i) = 1.0;
            vb[i] = y[i];
        }
        svd.solve();

        printf("y=%f x + %f\n\n", a, b);
        printf("a,b = ");
        print(svd.x());
        printf("\n");
        lmath::VectorView& vx = svd.x();
        lmath::element_type error = 0.0;
        for(lmath::s32 i=0; i<NumSamples; ++i){
            lmath::element_type e = vx[0]*x[i] + vx[1] - y[i];
            error += e*e;
        }
        printf("error = %f\n", error);

    }
    return 0;
}
