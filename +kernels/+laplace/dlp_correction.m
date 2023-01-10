function D_corr = dlp_correction(source, side, varargin)

% Parse optional inputs
p = inputParser;
addRequired(p, 'source', @(s) isa(s,'Boundary'));
addRequired(p, 'side', @(s) any(strcmp(s,{'i','e'})));
parse(p, source, side, varargin{:});

% Closeness factor
closepan = 1.2;

n = source.N;
np = source.np;
z = cell2mat(source.z);
kd = KDTree([real(z) imag(z)]);
D_corr = spalloc(n*np, n*np, 3*n*np);

%nf = 2*n;
nf = 32;
%nf = 16;
sourcef = resample(source, nf);
[x, ~, v] = legpts(n);
xf = legpts(nf);
%L = util.interpmat(x, xf);
L = barymat(xf, x, v);

for k = 1:np

    % Package source panel into structs for close evaluation routine
    % Original panel
    pa = struct();
    pa.x  = source.z{k};
    pa.xp = source.dz{k};
    pa.nx = source.normal{k}(:,1) + 1i*source.normal{k}(:,2);
    pa.sp = source.speed{k};
    pa.w  = source.w{k};
    pa.wxp = source.cw{k};

    % Upsampled panel
    paf = struct();
    paf.x  = sourcef.z{k};
    paf.xp = sourcef.dz{k};
    paf.nx = sourcef.normal{k}(:,1) + 1i*sourcef.normal{k}(:,2);
    paf.sp = sourcef.speed{k};
    paf.w  = sourcef.w{k};
    paf.wxp = sourcef.cw{k};

    % Find nearby panels using kd-tree
    a = source.zbreaks(k);
    b = source.zbreaks(k+1);
    center = (a+b)/2;
    center = [real(center) imag(center)];
    radius = closepan * abs(b-a);
    closeIdx = kd.ball(center, radius);

    % Don't correct self
    selfIdx = (1:n)+(k-1)*n;
    closeIdx = setdiff(closeIdx, selfIdx);

    % For each nearby panel, add the close evaluation matrix and subtract 
    % the direct evaluation matrix. This yields the correction matrix.
    if ( ~isempty(closeIdx) )
        t = struct();
        t.x = z(closeIdx);
        
        [~, idx] = min(abs(paf.x - t.x(:).'));
        v = t.x - paf.x(idx);
        w = paf.nx(idx);
        sgn = -0.5*(conj(v).*w + v.*conj(w)); % dot(v,w)
        D_close_i = LapDLP_closepanel(t, paf, a, b, 'i');
        D_close_e = LapDLP_closepanel(t, paf, a, b, 'e');
        D_close = D_close_i;
        D_close(sgn<0, :) = D_close_e(sgn<0, :);

        D_close  = LapDLP_closepanel(t, paf, a, b, side);
        D_direct = LapDLP(t, pa);
        D_close(isinf(D_direct)) = 0;
        D_direct(isinf(D_direct)) = 0;
        D_corr(closeIdx,selfIdx) = D_corr(closeIdx,selfIdx) + D_close*L - D_direct;
        %D_corr(closeIdx,selfIdx) = D_corr(closeIdx,selfIdx) + D_close*L;
    end

end

D_corr(1:size(D_corr,1)+1:end) = 0;

end
