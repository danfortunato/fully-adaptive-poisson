function u = SLP_eval_close_panel(source, target, sigma, side)

%u = SLP_eval_close_panel_ludvig(source, target, sigma, side);
%return

u = FMM_eval(source, target, 'charge', sigma);

n = source.N;
closepan = 1.2;
Z = target(:,1) + target(:,2)*1i;
for k = 1:source.np
    a = source.zbreaks(k);
    b = source.zbreaks(k+1);
    closeIdx = abs(Z - (a+b)/2) < closepan * abs(b-a);

    pa = struct();
    pa.x  = source.z{k};
    pa.xp = source.dz{k};
    pa.nx = source.normal{k}(:,1) + 1i*source.normal{k}(:,2);
    pa.sp = source.speed{k};
    pa.w  = source.w{k};
    pa.wxp = source.cw{k};

    if ( any(closeIdx) )
        t = struct();
        t.x = Z(closeIdx);
        S_close = LapSLP_closepanel(t, pa, a, b, side);
        S_direct = LapSLP(t, pa);
        den = sigma((1:n)+(k-1)*n);
        u(closeIdx) = u(closeIdx) + (S_close - S_direct) * den;
    end
end

end
