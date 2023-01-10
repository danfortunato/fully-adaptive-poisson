function u = DLP_eval_close_panel(source, target, sigma, side)

u = FMM_eval(source, target, 'dipstr', sigma);
%u=0;

n = source.N;
closepan = 1.2;
R = 5;

% 8th  order --> R = 8
% 12th order --> R = 4
% 16th order --> R = 3

kd = KDTree(target);

%nf = 32;
nf = 32;
%nf = 2*n;
%nf = 32;
sourcef = resample(source, nf);
[x, ~, v] = legpts(n);
xf = legpts(nf);
%L = util.interpmat(x, xf);
L = barymat(xf, x, v);

Z = target(:,1) + target(:,2)*1i;
%[~, idx] = sort(arclength(source, 1:source.np));
%idx = flip(idx);
for k = 1:source.np
%for k = idx(:).'

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
    %closeIdx = abs(Z - (a+b)/2) < closepan * abs(b-a);
    center = (a+b)/2;
    center = [real(center) imag(center)];
    radius = closepan * abs(b-a);
    closeIdx = kd.ball(center, radius);

%     %R = closepan * arclength(source, k);
%     closeIdxFiltered = [];
%     for j = 1:length(closeIdx)
%         tj = panel_preimage(pa.x, x, Z(closeIdx(j)));
%         %[bernstein_rad(tj) R]
%         if ( bernstein_rad(tj) < R )
%             closeIdxFiltered = [closeIdxFiltered; closeIdx(j)];
%         end
%     end
%     closeIdx = closeIdxFiltered;

    if ( any(closeIdx) )
    %if ( any(closeIdx) && (k == 200 || k == 201) )
    %if ( false )
        t = struct();
        t.x = Z(closeIdx);
        [~, idx] = min(abs(paf.x - t.x(:).'));
        v = t.x - paf.x(idx);
        w = paf.nx(idx);
        sgn = -0.5*(conj(v).*w + v.*conj(w)); % dot(v,w)
        ti = struct();
        te = struct();
        ti.x = t.x(sgn > 0);
        te.x = t.x(sgn < 0);
        %xa = t.x - a;
        %xb = t.x - b;
        %sgn = sign(imag(0.5*(conj(xa)*xb - xa*conj(xb))));
        %D_close = LapDLP_closepanel(t, paf, a, b, newside);
%         D_close_i = LapDLP_closepanel(t, paf, a, b, 'i');
%         D_close_e = LapDLP_closepanel(t, paf, a, b, 'e');
%         D_close = D_close_i;
%         D_close(sgn<0, :) = D_close_e(sgn<0, :);
        D_close = LapDLP_closepanel(t, paf, a, b, side);
        D_direct = LapDLP(t, pa);
        idx  = ( abs(t.x - pa.x.' ) < 1e-14 );
        idxf = ( abs(t.x - paf.x.') < 1e-14 );
        D_close(idxf) = 0;
        D_direct(idx) = 0;
        den = sigma((1:n)+(k-1)*n);
        u(closeIdx) = u(closeIdx) + (D_close*L - D_direct) * den;

%         z = source.z{k};
%         za = source.zbreaks(k);
%         zb = source.zbreaks(k+1);
%         dz = source.dz{k};
%         u_close = close_eval_dlp(z, s, dz, za, zb, den, Z(closeIdx));
%         u(closeIdx) = u(closeIdx) + u_close - D_direct*den;
    end
end

end

function rho = bernstein_rad(t)
    rho1 = abs(t + sqrt(t.^2-1));
    rho2 = abs(t - sqrt(t.^2-1));
    rho = max(rho1, rho2);
end

function rho = bernstein_rad2(z) 
    rho = abs(z + sqrt(z - 1).*sqrt(z+1));
end

function [t0, c, cp] = panel_preimage(z,t,z0)
% PANEL_PREIMAGE   solve for preimage of complex point under panel map.
%
% t0 = panel_preimage(z,t,z0) returns t0 preimage such that z0 = Z(t0) where
%  Z : C -> C is the analytic map that took the set of standard nodes t to
%  the set of nodes z.
%
% [t0 c cp] = panel_preimage(z,t,z0) also returns the coeffs (in polyval order)
%  of the approximation to the map Z and its complex derivative respectively.
%  c and cp do not depend on z0.

% to do: vectorize over z0, t0

% solve monomial rep of map via Vandermonde (indep of z0)
t = t(:); z = z(:); p = numel(t);
V = ones(p,p);
for k=2:p
  V(:,k) = t.*V(:,k-1);
end
c = V\z;                   % monomial coeffs
cp = c(2:end).*(1:p-1)';   % monomial coeffs of deriv of map
c=flipud(c); cp=flipud(cp);   % order max power to const term, for polyval

% Newton to solve for t0 in Z(t0)-z0 = 0.  todo: check vectorizes over z0, t0
maxit = 20;
zcen = (z(p)+z(1))/2; zsc = (z(p)-z(1))/(t(p)-t(1));
t0 = (z0-zcen)/zsc;    % initial guess
for i=1:maxit
  t0old = t0;
  t0 = t0 - (polyval(c,t0) - z0) ./ polyval(cp,t0);
  if max(abs(t0-t0old)) < 1e-15             % not rel since on [-1,1] scale
    break;
  end
end

end
