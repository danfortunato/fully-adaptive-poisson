#include <cmath>
#include "omp.h"
#include "sctl.hpp"

double norm(double x, double y)
{
    return std::sqrt(x*x + y*y);
}

template<int N>
double bary_vec(double x, const double fk[N], const double xk[N], const double vk[N])
{
    sctl::Vec<double,sctl::DefaultVecLen<double>()> x_vec(x);
    auto val_vec = sctl::Vec<double,sctl::DefaultVecLen<double>()>::Zero();
    auto sum_vec = sctl::Vec<double,sctl::DefaultVecLen<double>()>::Zero();
    for (int i=0; i<N; i+=sctl::DefaultVecLen<double>()) {
        auto xk_vec = sctl::Vec<double,sctl::DefaultVecLen<double>()>::LoadAligned(&xk[i]);
        auto vk_vec = sctl::Vec<double,sctl::DefaultVecLen<double>()>::LoadAligned(&vk[i]);
        auto fk_vec = sctl::Vec<double,sctl::DefaultVecLen<double>()>::LoadAligned(&fk[i]);
        sctl::Vec<double,sctl::DefaultVecLen<double>()> x_xk = x_vec - xk_vec;
        sctl::Vec<double,sctl::DefaultVecLen<double>()> xx = vk_vec * x_xk * sctl::approx_rsqrt<10,sctl::Vec<double,sctl::DefaultVecLen<double>()>>(x_xk*x_xk*x_xk*x_xk);//vk_vec / (x_vec - xk_vec);
        val_vec += xx * fk_vec;
        sum_vec += xx;
    }

    // Sum the partial sums
    double val = 0;
    double sum = 0;
    for (int i=0; i<sctl::DefaultVecLen<double>(); ++i) {
        val += val_vec[i];
        sum += sum_vec[i];
    }
    val /= sum;

    // Clean up NaN values
    if (std::isnan(val)) {
        for (int i=0; i<N; ++i) {
            if (x == xk[i]) {
                return fk[i];
            }
        }
    }

    return val;
}

template<int N>
double bary(double x, const double fk[N], const double xk[N], const double vk[N])
{
    double val = 0;
    double sum = 0;
    for (int i=0; i<N; ++i) {
        double xx = vk[i] / (x - xk[i]);
        val += xx * fk[i];
        sum += xx;
    }
    val /= sum;

    // Clean up NaN values
    if (std::isnan(val)) {
        for (int i=0; i<N; ++i) {
            if (x == xk[i]) {
                return fk[i];
            }
        }
    }

    return val;
}

template<int N>
void getElemCoordinates(int npts, int nthreads,
                        const double *x,     const double *y,
                        const double *gamx,  const double *gamy,
                        const double *gam1x, const double *gam1y,
                        const double *nx,    const double *ny,
                        const double *xcheb, const double *vcheb,
                        double *t, double *r, int *found)
{
    double tol = 1e-14;
    double perpx = -ny[0];
    double perpy =  nx[0];

    double det[N];
    for (int i=0; i<N; ++i) {
        det[i] = 1. / (perpy*nx[i] - perpx*ny[i]);
    }

    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int i=0; i<npts; ++i) {
        double qx = x[i];
        double qy = y[i];
        double tk[N];
        found[i] = 0;
        int sign_old = 0;

        for (int j=0; j<N; ++j) {
            double rhsx = gam1x[j] - qx;
            double rhsy = gam1y[j] - qy;
            double t_intersect = det[j] * (-perpy*rhsx + perpx*rhsy);
            double xk = gam1x[j] + nx[j]*t_intersect;
            double yk = gam1y[j] + ny[j]*t_intersect;
            double val = (qx-xk)*(qx-perpx) + (qy-yk)*(qy-perpy);
            val = (std::abs(val) < tol) ? 0 : val;
            int sign = (val > 0) - (val < 0);
            if (j>0 && sign != sign_old) {
                // The point is inside this element if the sign of intersection
                // has changed
                found[i] = 1;
            }
            sign_old = sign;
            tk[j] = sign * norm(qx-xk, qy-yk);
        }

        // Compute tangential coordinate
        t[i] = bary<N>(0, xcheb, tk, vcheb);
        t[i] = bary<N>(0, xcheb, tk, vcheb);

        // Evaluate tangential coordinate at the top and bottom of the element
        double gam_qx  = bary<N>(t[i], gamx,  xcheb, vcheb);
        double gam_qy  = bary<N>(t[i], gamy,  xcheb, vcheb);
        double gam1_qx = bary<N>(t[i], gam1x, xcheb, vcheb);
        double gam1_qy = bary<N>(t[i], gam1y, xcheb, vcheb);
        gam_qx  = bary<N>(t[i], gamx,  xcheb, vcheb);
        gam_qy  = bary<N>(t[i], gamy,  xcheb, vcheb);
        gam1_qx = bary<N>(t[i], gam1x, xcheb, vcheb);
        gam1_qy = bary<N>(t[i], gam1y, xcheb, vcheb);

        // Compute radial coordinate and map to [0,1]
        r[i] = 2*norm(qx-gam_qx, qy-gam_qy)/norm(gam1_qx-gam_qx, gam1_qy-gam_qy) - 1;
    }
}

// MWrap can't handle templates so we have to make a function for each N we want
template<typename... Args> void getElemCoordinates4(Args... args) { getElemCoordinates<4>(args...); }
template<typename... Args> void getElemCoordinates5(Args... args) { getElemCoordinates<5>(args...); }
template<typename... Args> void getElemCoordinates6(Args... args) { getElemCoordinates<6>(args...); }
template<typename... Args> void getElemCoordinates7(Args... args) { getElemCoordinates<7>(args...); }
template<typename... Args> void getElemCoordinates8(Args... args) { getElemCoordinates<8>(args...); }
template<typename... Args> void getElemCoordinates9(Args... args) { getElemCoordinates<9>(args...); }
template<typename... Args> void getElemCoordinates10(Args... args) { getElemCoordinates<10>(args...); }
template<typename... Args> void getElemCoordinates11(Args... args) { getElemCoordinates<11>(args...); }
template<typename... Args> void getElemCoordinates12(Args... args) { getElemCoordinates<12>(args...); }
template<typename... Args> void getElemCoordinates13(Args... args) { getElemCoordinates<13>(args...); }
template<typename... Args> void getElemCoordinates14(Args... args) { getElemCoordinates<14>(args...); }
template<typename... Args> void getElemCoordinates15(Args... args) { getElemCoordinates<15>(args...); }
template<typename... Args> void getElemCoordinates16(Args... args) { getElemCoordinates<16>(args...); }
template<typename... Args> void getElemCoordinates17(Args... args) { getElemCoordinates<17>(args...); }
template<typename... Args> void getElemCoordinates18(Args... args) { getElemCoordinates<18>(args...); }
template<typename... Args> void getElemCoordinates19(Args... args) { getElemCoordinates<19>(args...); }
template<typename... Args> void getElemCoordinates20(Args... args) { getElemCoordinates<20>(args...); }
template<typename... Args> void getElemCoordinates21(Args... args) { getElemCoordinates<21>(args...); }
template<typename... Args> void getElemCoordinates22(Args... args) { getElemCoordinates<22>(args...); }
template<typename... Args> void getElemCoordinates23(Args... args) { getElemCoordinates<23>(args...); }
template<typename... Args> void getElemCoordinates24(Args... args) { getElemCoordinates<24>(args...); }
template<typename... Args> void getElemCoordinates25(Args... args) { getElemCoordinates<25>(args...); }
template<typename... Args> void getElemCoordinates26(Args... args) { getElemCoordinates<26>(args...); }
template<typename... Args> void getElemCoordinates27(Args... args) { getElemCoordinates<27>(args...); }
template<typename... Args> void getElemCoordinates28(Args... args) { getElemCoordinates<28>(args...); }
template<typename... Args> void getElemCoordinates29(Args... args) { getElemCoordinates<29>(args...); }
template<typename... Args> void getElemCoordinates30(Args... args) { getElemCoordinates<30>(args...); }
template<typename... Args> void getElemCoordinates31(Args... args) { getElemCoordinates<31>(args...); }
template<typename... Args> void getElemCoordinates32(Args... args) { getElemCoordinates<32>(args...); }
