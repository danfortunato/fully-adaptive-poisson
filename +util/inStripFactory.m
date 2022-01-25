function inStrip = inStripFactory(Gamma, Gamma1)

    sigma = ones(Gamma.N*Gamma.np, 1);

    function inStripFunc = inStripFunc(xy)
        outGamma = kernels.laplace.dlp(Gamma,  'target', xy, 'density', sigma, 'side', 'e');
        inGamma1 = kernels.laplace.dlp(Gamma1, 'target', xy, 'density', sigma, 'side', 'i');
        inStripFunc = abs(outGamma) > 0.5 & abs(inGamma1+1) > 0.5;
    end

inStrip = @inStripFunc;

end
