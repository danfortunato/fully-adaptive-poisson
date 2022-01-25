function u = DLP_eval_close_panel(source, target, sigma, side)

%u = DLP_eval_close_panel_ludvig(source, target, sigma, side);
%u = DLP_eval_close_panel_old(source, target, sigma, side);
%return

u = FMM_eval(source, target, 'dipstr', sigma);

s = legpts(source.N);

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
        D_close = LapDLP_closepanel(t, pa, a, b, side);
        D_direct = LapDLP(t, pa);
        den = sigma((1:n)+(k-1)*n);
        u(closeIdx) = u(closeIdx) + (D_close - D_direct) * den;

%         z = source.z{k};
%         za = source.zbreaks(k);
%         zb = source.zbreaks(k+1);
%         dz = source.dz{k};
%         u_close = close_eval_dlp(z, s, dz, za, zb, den, Z(closeIdx));
%         u(closeIdx) = u(closeIdx) + u_close - D_direct*den;
    end
end

end
