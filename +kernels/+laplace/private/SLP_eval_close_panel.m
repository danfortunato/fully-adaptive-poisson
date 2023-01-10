function u = SLP_eval_close_panel(source, target, sigma, side)

u = FMM_eval(source, target, 'charge', sigma);

n = source.N;
closepan = 1;

% 8th  order --> R = 8
% 12th order --> R = 4
% 16th order --> R = 3

kd = KDTree(target);

%nf = 32;
%nf = 16;
%nf = 2*n;
nf = 32;
sourcef = resample(source, nf);
[x, ~, v] = legpts(n);
xf = legpts(nf);
%L = util.interpmat(x, xf);
L = barymat(xf, x, v);

Z = target(:,1) + target(:,2)*1i;
for k = 1:source.np

    pa = struct();
    pa.x  = source.z{k};
    pa.xp = source.dz{k};
    pa.nx = source.normal{k}(:,1) + 1i*source.normal{k}(:,2);
    pa.sp = source.speed{k};
    pa.w  = source.w{k};
    pa.wxp = source.cw{k};

    paf = struct();
    paf.x  = sourcef.z{k};
    paf.xp = sourcef.dz{k};
    paf.nx = sourcef.normal{k}(:,1) + 1i*sourcef.normal{k}(:,2);
    paf.sp = sourcef.speed{k};
    paf.w  = sourcef.w{k};
    paf.wxp = sourcef.cw{k};

    a = source.zbreaks(k);
    b = source.zbreaks(k+1);
    center = (a+b)/2;
    center = [real(center) imag(center)];
    radius = closepan * abs(b-a);
    closeIdx = kd.ball(center, radius);

    if ( any(closeIdx) )
        t = struct();
        t.x = Z(closeIdx);
        S_close = LapSLP_closepanel(t, paf, a, b, side);
        S_direct = LapSLP(t, pa);
        idx  = ( abs(t.x - pa.x.' ) < 1e-14 );
        idxf = ( abs(t.x - paf.x.') < 1e-14 );
        S_close(idxf) = 0;
        S_direct(idx) = 0;
        den = sigma((1:n)+(k-1)*n);
        u(closeIdx) = u(closeIdx) + (S_close*L - S_direct) * den;
    end

end

end
