function C = levelrestrict(B, varargin)

% Parse optional inputs
p = inputParser;
addRequired(p, 'B', @(s) isa(s,'Boundary'));
addParameter(p, 'fac', 2.1, @isnumeric);
addParameter(p, 'maxsize', Inf, @isnumeric);
parse(p, B, varargin{:});
fac     = p.Results.fac;
maxsize = p.Results.maxsize;

kind = 'arclength';
npmax = 10000;
maxiter = 1000;
splititer = 20;
np = B.np;

% Make an arc length parametrization so that panels are divided into two
% pieces of equal arc length
%B = reparametrize(B);
[x01, w01, v] = chebpts(B.N, [0 1]);
w01 = w01(:).';
x = cell(np, 1);
w = cell(np, 1);
speed = cell(np, 1);
ints = cell(np, 1);
I = intmat(B.N, 1, [0 1]);
for k = 1:np
    x{k} = (B.breaks(k+1)-B.breaks(k))*x01 + B.breaks(k);
    w{k} = (B.breaks(k+1)-B.breaks(k))*w01;
    speed{k} = legvals2chebvals(B.speed{k});
    %ints{k} = intmat(B.N, 1, B.breaks(k:k+1)) * speed{k};
    ints{k} = diff(B.breaks(k:k+1)) * (I * speed{k});
end

% [s, sbreaks] = arclengthparam(B);
% for k = 1:np
%     %x{k} = (sbreaks(k+1)-sbreaks(k))*x01 + sbreaks(k);
%     %w{k} = (sbreaks(k+1)-sbreaks(k))*w01;
%     x{k} = (B.breaks(k+1)-B.breaks(k))*x01 + B.breaks(k);
%     w{k} = (B.breaks(k+1)-B.breaks(k))*w01;
%     ints{k} = legvals2chebvals(s(:,k)) - sbreaks(k);
% end

ab = zeros(2,npmax);
ab(1,1:np) = B.breaks(1:end-1);
ab(2,1:np) = B.breaks(2:end);
%ab(1,1:np) = sbreaks(1:end-1);
%ab(2,1:np) = sbreaks(2:end);

adjs = zeros(2,npmax);
adjs(1,1:np) = [np 1:np-1];
adjs(2,1:np) = [2:np 1];

parents = zeros(npmax, 1);
parents(1:np) = 1:np;

maxself = 0;

iter = 0;
done = false;
while ( ~done && iter < maxiter )
    %[np done iter]

    done = true;
    npold = np;
    %P = randperm(npold);
    P = 1:npold;
    [~, idx] = sort(diff(ab(:,1:npold)));
    P = P(idx);
    for k = P(:).'
        kl = adjs(1,k);
        kr = adjs(2,k);

        a = ab(1,k);
        b = ab(2,k);

        if ( strcmpi(kind, 'arclength') )
            self  = diff(bary(ab(:,k),  ints{parents(k)},  x{parents(k)},  v));
            left  = diff(bary(ab(:,kl), ints{parents(kl)}, x{parents(kl)}, v));
            right = diff(bary(ab(:,kr), ints{parents(kr)}, x{parents(kr)}, v));
        else
            self  = b-a;
            left  = ab(2,kl)-ab(1,kl);
            right = ab(2,kr)-ab(1,kr);
        end

        maxself = max(self, maxself);

        if ( self > fac*left || self > fac*right || self > maxsize )
            %[self left right]

            if ( np == npmax )
                error('Too many panels.');
            end

            done = false;

            % Do a few iterations of binary search to find a crude estimate
            % of the arc-length midpoint:

            alen = bary(a, ints{parents(k)}, x{parents(k)}, v);
            blen = bary(b, ints{parents(k)}, x{parents(k)}, v);
            targetlen = (blen-alen)/2;
            lower = a;
            upper = b;
            for l = 1:splititer
                mid = (upper+lower)/2;
                midlen = bary(mid, ints{parents(k)}, x{parents(k)}, v) - alen;
                if ( midlen < targetlen )
                    lower = mid;
                else
                    upper = mid;
                end
            end
            %mid = (a+b)/2;

            if ( ~(a < mid || mid < b) )
                keyboard
            end

            knew = np+1;
            adjs(2,k) = knew;
            adjs(:,knew) = [k kr];
            adjs(1,kr) = knew;
            parents(knew) = parents(k);

            ab(:,k)    = [a mid];
            ab(:,knew) = [mid b];

            np = np+1;
        end

    end

    iter = iter+1;

end

%ab = ab(:,1:np);2
%adjs = adjs(:,1:np);

breaks = [sort(ab(1,1:np),2), B.breaks(end)];
%breaks = [sort(ab(1,1:np),2), sbreaks(end)];

[s01, w01] = util.gauss(B.N, 0, 1);
np = numel(breaks)-1;
% s = cell(np,1);
% w = cell(np,1);
% z = cell(np,1);
% for k = 1:np
%     s{k} = (breaks(k+1)-breaks(k))*s01 + breaks(k);
%     w{k} = (breaks(k+1)-breaks(k))*w01;
%     z{k} = B.f(s{k});
% end
% 
% quadinfo = [];
% quadinfo.s = s;
% quadinfo.w = w;
% quadinfo.breaks = breaks;
% 
% C = Boundary(z, quadinfo);

C = Boundary();
C.f = B.f;
C.df = B.df;
C.dff = B.dff;
C.bend = B.bend;
C.levelset = B.levelset;
C.rule = B.rule;
C.N = B.N;

C.np = np;
C.breaks = breaks;
C.s         = cell(C.np,1);
C.z         = cell(C.np,1);
C.dz        = cell(C.np,1);
C.dzz       = cell(C.np,1);
C.x         = cell(C.np,1);
C.dx        = cell(C.np,1);
C.dxx       = cell(C.np,1);
C.w         = cell(C.np,1);
C.cw        = cell(C.np,1);
C.speed     = cell(C.np,1);
C.accel     = cell(C.np,1);
C.normal    = cell(C.np,1);
C.curvature = cell(C.np,1);

j = 1;
[~, ~, D] = util.gauss(C.N);
for k = 1:C.np
    C.s{k} = (C.breaks(k+1)-C.breaks(k))*s01 + breaks(k);

%     if ( ~(C.breaks(k) >= B.breaks(j) && C.breaks(k+1) <= B.breaks(j+1)) )
%         j = j+1;
%     end
%     [~, ~, v01] = legpts(B.N, B.breaks(j:j+1));
%     C.z{k}   = bary(C.s{k}, B.z{j},   B.s{j}, v01);
%     C.dz{k}  = bary(C.s{k}, B.dz{j},  B.s{j}, v01);
%     C.dzz{k} = bary(C.s{k}, B.dzz{j}, B.s{j}, v01);
%     C.zbreaks(k:k+1) = bary(C.breaks(k:k+1), B.z{j}, B.s{j}, v01);
    %C.speed{k} = bary(C.s{k}, B.speed{j}, B.s{j}, v01);
    %C.accel{k} = bary(C.s{k}, B.accel{j}, B.s{j}, v01);
    %C.normal{k} = bary(C.s{k}, B.normal{j}, B.s{j}, v01);
    %C.curvature{k} = bary(C.s{k}, B.curvature{j}, B.s{j}, v01);

    C.z{k} = C.f(C.s{k});
    scl = 2 / (C.breaks(k+1) - C.breaks(k));
    if ( ~isempty(C.df) )
        C.dz{k}  = C.df(C.s{k});
    else
        C.dz{k} = scl * (D * C.z{k});
    end
    if ( ~isempty(C.dff) )
        C.dzz{k} = C.dff(C.s{k});
    else
        C.dzz{k} = scl * (D * C.dz{k});
    end

    C.x{k}   = [ real(C.z{k})   imag(C.z{k})   ];
    C.dx{k}  = [ real(C.dz{k})  imag(C.dz{k})  ];
    C.dxx{k} = [ real(C.dzz{k}) imag(C.dzz{k}) ];
    C.speed{k} = sqrt(sum(C.dx{k}.^2, 2));
    C.accel{k} = real(C.dz{k} ./ C.speed{k} .* conj(C.dzz{k}));
    C.normal{k} = [ C.dx{k}(:,2), -C.dx{k}(:,1) ] ./ C.speed{k};
    C.curvature{k} = -dot(C.dxx{k}, C.normal{k}, 2) ./ C.speed{k}.^2;
    w = (C.breaks(k+1)-C.breaks(k))*w01;
    C.w{k}  = w .* C.speed{k};
    C.cw{k} = w .* C.dz{k};
end

%[xleg, ~, vleg] = legpts(C.N);
%C.zbreaks = bary(-1, [C.z{:}], xleg, vleg);
%C.zbreaks(end+1) = C.zbreaks(1);
C.zbreaks = C.f(C.breaks);
                
nodes = cell2mat(C.x);
npts = size(nodes, 1);
C.polynodes = nodes;
C.polyedges = [1:npts ; 2:npts 1].';

end
