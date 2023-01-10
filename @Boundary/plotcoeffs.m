function plotcoeffs(dom)

n = dom.N;
for k = 1:length(dom)
    cfs = legvals2legcoeffs(dom.z{k});
    semilogy(0:1:n-1, abs(cfs), '.', 'MarkerSize', 20)
    hold on
end

end
