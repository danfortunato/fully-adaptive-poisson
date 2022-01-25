//#define EIGEN_STACK_ALLOCATION_LIMIT 10000000
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_USE_MKL_ALL
#define MKL_ALIGN 64

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <mkl.h>
#include "omp.h"

using Eigen::all;
using Eigen::seq;

template<int M = Eigen::Dynamic, int N = M> using Matrix  = Eigen::Matrix<double, M, N>;
template<int M = Eigen::Dynamic, int N = M> using Array   = Eigen::Array<double, M, N>;
template<int M = Eigen::Dynamic, int N = M> using Indices = Eigen::Array<int, M, N>;
template<typename T> using Map = Eigen::Map<T>;

template<int N, typename DerivedA, typename DerivedB>
Matrix<> kron(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedB>& B)
{
    Matrix<> C(N*N, N*N);
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            C.template block<N,N>(N*i,N*j) = A(i,j) * B;
        }
    }
    return C;
}

template<int N, typename Derived>
void normals(const Eigen::MatrixBase<Derived>& x,
             const Eigen::MatrixBase<Derived>& y,
             const Eigen::MatrixBase<Derived>& D,
             Array<N-2,2>& nl,
             Array<N-2,2>& nr,
             Array<N-2,2>& nd,
             Array<N-2,2>& nu)
{
    Array<N,1> dx, dy, scl;

    dx = D * x(all,0).reverse(); // reverse?
    dy = D * y(all,0).reverse();
    scl = (dx*dx + dy*dy).sqrt();
    dx = dx / scl;
    dy = dx / scl;
    nl << -dy.template segment<N-2>(1), dx.template segment<N-2>(1);
    //nl = -[dy, -dx] / (dx*dx + dy*dy).sqrt();
    //nl = nl(seq(1, n-2), all).eval(); // nl([1,n],:) = [];

    dx = D * x(all,N-1).reverse();
    dy = D * y(all,N-1).reverse();
    scl = (dx*dx + dy*dy).sqrt();
    dx = dx / scl;
    dy = dx / scl;
    nr << dy.template segment<N-2>(1), -dx.template segment<N-2>(1);

    dx = D*x(0,all).transpose();
    dy = D*y(0,all).transpose();
    scl = (dx*dx + dy*dy).sqrt();
    dx = dx / scl;
    dy = dx / scl;
    nd << -dy.template segment<N-2>(1), dx.template segment<N-2>(1);

    dx = D*x(N-1,all).transpose();
    dy = D*y(N-1,all).transpose();
    scl = (dx*dx + dy*dy).sqrt();
    dx = dx / scl;
    dy = dx / scl;
    nu << dy.template segment<N-2>(1), -dx.template segment<N-2>(1);
}

template<int N, typename Derived>
void curvedDiffmat(const Eigen::MatrixBase<Derived>& x,
                   const Eigen::MatrixBase<Derived>& y,
                   const Eigen::MatrixBase<Derived>& D, 
                   Matrix<>& Dx,
                   Matrix<>& Dy)
{
    // Compute the Jacobian:
    Array<N> dxdt = (x * D.transpose());//.array();
    Array<N> dydt = (y * D.transpose()).array();
    Array<N> dxdr = (D * x).array();
    Array<N> dydr = (D * y).array();
    Array<N> J = -(dxdr*dydt - dxdt*dydr);

    // Compute d[t,r]/d[x,y] as the inverse of the Jacobian:
    Array<N> dtdx = -dydr / J;
    Array<N> dtdy =  dxdr / J;
    Array<N> drdx =  dydt / J;
    Array<N> drdy = -dxdt / J;

    // Construct curved differentiation matrices:
    auto I = Matrix<N>::Identity();
    Matrix<> Dt2d = kron<N>(D, I);
    Matrix<> Dr2d = kron<N>(I, D);
    Dx = Dt2d.array().colwise()*dtdx.reshaped() + Dr2d.array().colwise()*drdx.reshaped();
    Dy = Dt2d.array().colwise()*dtdy.reshaped() + Dr2d.array().colwise()*drdy.reshaped();
}

template<int N>
void constructOperators(int numPatches,
                        int nthreads,
                        const double* x_data,
                        const double* y_data,
                        const double* rhs_data,
                        const double* D_data,
                        const double* B_data,
                        double* S_data,
                        double* D2N_data,
                        double* Aii_data)
{
    const int numIntPts = (N-2)*(N-2) + 4;
    const int numBdyPts = 4*(N-2);

    // Chebyshev differentiation matrix
    Map<const Matrix<N>> D(D_data);

    // Corner interpolation matrix
    Map<const Matrix<2,N-2>> B(B_data);

    // Interior and exterior indices
    Array<N> I = Array<N>::Zero();
    I.template block<N-2,N-2>(1,1) = 1;
    I({0,N-1},{0,N-1}) = 1;
    Indices<Eigen::Dynamic,1> ii(numIntPts);
    Indices<Eigen::Dynamic,1> ee(numBdyPts);
    int ki=0, ke=0;
    for (int ij=0; ij<N*N; ij++) {
        if (I(ij) == 1) {
            ii(ki) = ij;
            ki++;
        } else {
            ee(ke) = ij;
            ke++;
        }
    }

    // Indices for each side of the patch
    auto leftIdx  = seq(0,       1*(N-2)-1);
    auto rightIdx = seq(N-2,     2*(N-2)-1);
    auto downIdx  = seq(2*(N-2), 3*(N-2)-1);
    auto upIdx    = seq(3*(N-2), 4*(N-2)-1);
    int ibc = 3*(N-2)+1;

    // Scattering indices for the operators
    Indices<numBdyPts, 1> ss;
    int k = 0;
    for (int i=0;     i<N-2;     i++,  k++) ss(k) = i;
    for (int i=ibc-1; i<4*(N-2); i++,  k++) ss(k) = i;
    for (int i=N-2;   i<ibc-1;   i+=2, k++) ss(k) = i;
    for (int i=N-1;   i<ibc-1;   i+=2, k++) ss(k) = i;

    #pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int i=0; i<numPatches; i++) {

        // Get the data for this patch
        Map<const Matrix<N>> x(&x_data[N*N*i]), y(&y_data[N*N*i]);
        Map<const Matrix<numIntPts,1>> rhs(&rhs_data[numIntPts*i]);
        Map<Matrix<N*N,numBdyPts+1>> S(&S_data[N*N*(numBdyPts+1)*i]);
        Map<Matrix<numBdyPts,numBdyPts+1>> D2N(&D2N_data[numBdyPts*(numBdyPts+1)*i]);
        Map<Matrix<numIntPts>> Aii(&Aii_data[numIntPts*numIntPts*i]);

        // Construct curved differentiation matrices
        Matrix<> Dx(N*N, N*N), Dy(N*N, N*N);
        curvedDiffmat<N>(x, y, D, Dx, Dy);

        // Construct solution operator
        Matrix<> A = Dx*Dx + Dy*Dy;
        //Eigen::MatrixXd Aii2 = A(ii,ii).eval();
        Aii = A(ii,ii);
        Matrix<> F(numIntPts, numBdyPts+1);
        F << -A(ii,ee), rhs;
        //Matrix<> Sii = Aii.eval();//.llt().solve(F);
        Matrix<> Sii = Aii.llt().solve(F);

        // // Replace solution operator for corners with interp conditions:
        // Sii({0,1,numIntPts-2,numIntPts-1}, all).setZero();                  // S([1:2,end-1:end],:) = 0;
        // Sii({0,1}, seq(0,N-3)) = B;                                         // S(1:2,1:n-2) = B;
        // Sii({numIntPts-2,numIntPts-1}, seq(numBdyPts-N+2,numBdyPts-1)) = B; // S([end-1,end],end-n+2:end-1) = B;

        // // Append boundary points to solution operator:
        // Matrix tmpS(N*N, numBdyPts+1);
        // tmpS(ii, all) = Sii;
        // tmpS(ee, all).setIdentity();
        // Indices<numBdyPts+1, 1> idx;
        // idx << ss, numBdyPts;
        // S = tmpS(all, idx);

        // // Construct the D2N map:
        // //Matrix<numBdyPts, numBdyPts+1> D2N;
        // Array<N-2,2> nl, nr, nd, nu;
        // normals<N>(x, y, D, nl, nr, nd, nu);
        // Array<numBdyPts, numBdyPts+1> dx = Dx(ee, all).eval() * S;
        // Array<numBdyPts, numBdyPts+1> dy = Dy(ee, all).eval() * S;
        // dx = dx(ss, all);
        // dy = dy(ss, all);
        // D2N(leftIdx,  all) = dx(leftIdx,  all).colwise()*nl(all,0) + dy(leftIdx,  all).colwise()*nl(all,1);
        // D2N(rightIdx, all) = dx(rightIdx, all).colwise()*nr(all,0) + dy(rightIdx, all).colwise()*nr(all,1);
        // D2N(downIdx,  all) = dx(downIdx,  all).colwise()*nd(all,0) + dy(downIdx,  all).colwise()*nd(all,1);
        // D2N(upIdx,    all) = dx(upIdx,    all).colwise()*nu(all,0) + dy(upIdx,    all).colwise()*nu(all,1);
    }
}

// MWrap can't handle templates so we have to make a function for each N we want
template<typename... Args> void constructOperators4(Args... args) { constructOperators<4>(args...); }
template<typename... Args> void constructOperators5(Args... args) { constructOperators<5>(args...); }
template<typename... Args> void constructOperators6(Args... args) { constructOperators<6>(args...); }
template<typename... Args> void constructOperators7(Args... args) { constructOperators<7>(args...); }
template<typename... Args> void constructOperators8(Args... args) { constructOperators<8>(args...); }
template<typename... Args> void constructOperators9(Args... args) { constructOperators<9>(args...); }
template<typename... Args> void constructOperators10(Args... args) { constructOperators<10>(args...); }
template<typename... Args> void constructOperators11(Args... args) { constructOperators<11>(args...); }
template<typename... Args> void constructOperators12(Args... args) { constructOperators<12>(args...); }
template<typename... Args> void constructOperators13(Args... args) { constructOperators<13>(args...); }
template<typename... Args> void constructOperators14(Args... args) { constructOperators<14>(args...); }
template<typename... Args> void constructOperators15(Args... args) { constructOperators<15>(args...); }
template<typename... Args> void constructOperators16(Args... args) { constructOperators<16>(args...); }
template<typename... Args> void constructOperators17(Args... args) { constructOperators<17>(args...); }
template<typename... Args> void constructOperators18(Args... args) { constructOperators<18>(args...); }
template<typename... Args> void constructOperators19(Args... args) { constructOperators<19>(args...); }
template<typename... Args> void constructOperators20(Args... args) { constructOperators<20>(args...); }
template<typename... Args> void constructOperators21(Args... args) { constructOperators<21>(args...); }
template<typename... Args> void constructOperators22(Args... args) { constructOperators<22>(args...); }
template<typename... Args> void constructOperators23(Args... args) { constructOperators<23>(args...); }
template<typename... Args> void constructOperators24(Args... args) { constructOperators<24>(args...); }
template<typename... Args> void constructOperators25(Args... args) { constructOperators<25>(args...); }
template<typename... Args> void constructOperators26(Args... args) { constructOperators<26>(args...); }
template<typename... Args> void constructOperators27(Args... args) { constructOperators<27>(args...); }
template<typename... Args> void constructOperators28(Args... args) { constructOperators<28>(args...); }
template<typename... Args> void constructOperators29(Args... args) { constructOperators<29>(args...); }
template<typename... Args> void constructOperators30(Args... args) { constructOperators<30>(args...); }
template<typename... Args> void constructOperators31(Args... args) { constructOperators<31>(args...); }
template<typename... Args> void constructOperators32(Args... args) { constructOperators<32>(args...); }

