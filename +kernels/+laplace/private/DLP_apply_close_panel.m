function u = DLP_apply_close_panel(source, sigma, side)

target = cell2mat(source.x);
u = FMM_eval(source, target, 'dipstr', sigma);

n = source.N;
closepan = 1.2;
%closepan = 2.5;

nf = 2*n;
sourcef = resample(source, nf);
[x, ~, v] = legpts(n);
xf = legpts(nf);
%L = util.interpmat(x, xf);
L = barymat(xf, x, v);

Z_all = target(:,1) + target(:,2)*1i;
for k = 1:source.np
    
    ak = source.zbreaks(k);
    bk = source.zbreaks(k+1);
    midk = (ak+bk)/2;
    
    pa = struct();
    pa.x  = source.z{k};
    pa.xp = source.dz{k};
    pa.nx = source.normal{k}(:,1) + 1i*source.normal{k}(:,2);
    pa.sp = source.speed{k};
    pa.w  = source.w{k};
    pa.wxp = source.cw{k};
    
    den = sigma((1:n)+(k-1)*n);
    
    for j = 1:source.np
        if ( j == k ), continue, end
        aj = source.zbreaks(j);
        bj = source.zbreaks(j+1);
        midj = (aj+bj)/2;
        isclose = any(abs(midj - midk) < closepan * abs(bk-ak));

        if ( isclose )
            t = struct();
            t.x = source.z{j};
            D_close = LapDLP_closepanel(t, pa, ak, bk, side);
            D_direct = LapDLP(t, pa);
            jidx = (j-1)*n+(1:n);
            u(jidx) = u(jidx) + (D_close - D_direct) * den;
        end
    end

%     Z = Z_all;
%     Z((k-1)*n+(1:n)) = [];
%     a = source.zbreaks(k);
%     b = source.zbreaks(k+1);
%     closeIdx = abs(Z - (a+b)/2) < closepan * abs(b-a);
% 
%     pa = struct();
%     pa.x  = source.z{k};
%     pa.xp = source.dz{k};
%     pa.nx = source.normal{k}(:,1) + 1i*source.normal{k}(:,2);
%     pa.sp = source.speed{k};
%     pa.w  = source.w{k};
%     pa.wxp = source.cw{k};
% 
%     if ( any(closeIdx) )
%         t = struct();
%         t.x = Z(closeIdx);
%         D_close = LapDLP_closepanel(t, pa, a, b, side);
%         D_direct = LapDLP(t, pa);
%         den = sigma((1:n)+(k-1)*n);
%         u(closeIdx) = u(closeIdx) + (D_close - D_direct) * den;
%     end
end

end
