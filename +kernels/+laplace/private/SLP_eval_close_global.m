function u = SLP_eval_close_global(source, target, sigma, side)

x = cell2mat(source.x);
z = cell2mat(source.z);
s = cell2mat(source.s);
cw = cell2mat(source.cw);
w = cell2mat(source.w);
a = mean(z);

% Step 1: eval v^+ or v^- = cmplx SLP(tau):
vb = CSLPselfmatrix(source, side) * sigma;
%vb = SLP_nystrom(source) * sigma;
%vb = FMM_eval(source, x, 'charge', sigma);
if ( side == 'e' )
    sawlog = s/1i + log(a - z);  % sawtooth with jump cancelled by log
    for i=1:numel(sawlog)-1
        p = imag(sawlog(i+1)-sawlog(i)); % remove phase jumps
        sawlog(i+1) = sawlog(i+1) - 2i*pi*round(p/(2*pi));
    end
    totchgp = sum(w.*sigma)/(2*pi); % total charge due to SLP, over 2pi
    vb = vb + sawlog * totchgp;     % is now r
    %cw = 1i*s.nx.*s.w;  % complex speed weights for native quadr...
    %vinf = sum(bsxfun(@times, vb, cw./(s.x-s.a)),1) / (2i*pi); % interior Cauchy gets v_infty
    vinf = sum(1./(z-a).*vb.*cw) / (2i*pi); % interior Cauchy gets v_infty
    %vb = vb - ones(size(vb,1),1)*vinf;  % kill off v_infty so that v is in exterior Hardy space
    vb = vb - vinf;
end

% Step 2: compensated close-evaluation of v & v', followed by take Re:
v = Cauchy_close_global(source, target, vb, side); % does Sec. 3 of [lsc2d]
%v = v*2;
u = real(v);

if ( side == 'e' )
    zt = target(:,1)+target(:,2)*1i;
    % add real part of log and of v_infty back in...
    %u = u - log(abs(a - zt))*totchgp + ones(M,1)*real(vinf);   % Re v_infty = 0 anyway
    u = u - log(abs(a - zt))*totchgp + real(vinf);   % Re v_infty = 0 anyway
end

end
