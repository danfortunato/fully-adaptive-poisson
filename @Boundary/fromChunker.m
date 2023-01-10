function C = fromChunker(chnkr)

chnkr = sort(chnkr);

x = cell(chnkr.nch, 1);
for k = 1:chnkr.nch
    x{k} = chnkr.r(:,:,k).';
end

C = Boundary(x);
C = reparametrize(C);

end
