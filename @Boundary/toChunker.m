function chnkr = toChunker(Gamma)

pref = [];
pref.k = Gamma.N;
chnkr = chunker(pref);
chnkr = chnkr.addchunk(Gamma.np);

for k = 1:Gamma.np
    chnkr.r(1,:,k) = Gamma.x{k}(:,1);
    chnkr.r(2,:,k) = Gamma.x{k}(:,2);
    chnkr.d(1,:,k) = Gamma.dx{k}(:,1);
    chnkr.d(2,:,k) = Gamma.dx{k}(:,2);
    chnkr.d2(1,:,k) = Gamma.dxx{k}(:,1);
    chnkr.d2(2,:,k) = Gamma.dxx{k}(:,2);
    chnkr.h(k) = diff(Gamma.breaks([k k+1])) / 2;
    chnkr.adj(:,k) = [k-1 k+1];
end
chnkr.adj(1,1) = Gamma.np;
chnkr.adj(2,end) = 1;
chnkr.breaks = Gamma.breaks(1:end-1);

end
