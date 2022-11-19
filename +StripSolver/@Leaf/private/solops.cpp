#include "solops.h"

int main()
{
    const int N = 32;
    int numPatches = 10;
    int nthreads = 1;
    std::vector<double> x(numPatches*N*N, 0);
    std::vector<double> y(numPatches*N*N, 0);
    std::vector<double> rhs(numPatches*N*N, 0);
    std::vector<double> D(N*N, 0);
    std::vector<double> B(2*(N-2), 0);

    const int numIntPts = (N-2)*(N-2) + 4;
    const int numBdyPts = 4*(N-2);
    std::vector<double> S(numPatches*N*N*(numBdyPts+1), 0);
    std::vector<double> D2N(numPatches*numBdyPts*(numBdyPts+1), 0);
    std::vector<double> Aii(numPatches*numIntPts*numIntPts, 0);

    constructOperators<N>(numPatches, nthreads, x.data(), y.data(), rhs.data(), D.data(), B.data(),
        S.data(), D2N.data(), Aii.data());
}
