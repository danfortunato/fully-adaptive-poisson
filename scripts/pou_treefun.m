%% Define Gamma
n = 32;      % Number of nodes per panel
m = n;       % Number of radial Chebyshev nodes in strip
%psem = 2*n;
psem = n;

Gamma = Boundary.circle(n, 'quadrature', 'panel');
Gamma = refine(Gamma, 2);
%Gamma = Boundary.star(n, 'quadrature', 'panel');
%Gamma = refine(Gamma);

gamxy = cell2mat(Gamma.x);
gamx  = gamxy(:,1);
gamy  = gamxy(:,2);

scale = 1.2;
dom_global = boundingbox(Gamma, scale);

sol = @(x,y) -cos(x).*exp(sin(x)).*sin(y);
f   = @(x,y) (2*cos(x)+3*cos(x).*sin(x)-cos(x).^3).*exp(sin(x)).*sin(y);
g   = @(x,y) sol(x,y);

f = @(x,y) 1+0*x;

%% Define the fake interior boundary, Gamma'
tic
beta = 4;
bleed = 10;
width = 0.5;
%width = 1;
[x, y] = smoothStrip2(Gamma, n, beta, bleed, width);
xleg = chebvals2legvals(reshape(x, n, Gamma.np)); xleg = xleg(:);
yleg = chebvals2legvals(reshape(y, n, Gamma.np)); yleg = yleg(:);
z1 = mat2cell(xleg + 1i*yleg, repmat(n, Gamma.np, 1), 1);
%z1 = mat2cell(xleg + 1i*yleg + 0.15, repmat(n, Gamma.np, 1), 1);
Gamma1 = Boundary(z1);

% zz = cell2mat(Gamma.z);
% normals = cell2mat(Gamma.normal);
% th = 0.5;
% %th = 0;
% R = [cos(th) -sin(th);
%      sin(th)  cos(th)];
% normals = (R*normals.').';
% cn = normals(:,1) + 1i*normals(:,2);
% zz1 = zz - 0.1*cn;
% z1 = mat2cell(zz1, repmat(n, Gamma.np, 1), 1);
% Gamma1 = Boundary(z1);

% Build the strip grid
t = chebpts(psem, [0 1]);
xsem = chebpts(psem);
xx = cell(Gamma.np,1);
yy = cell(Gamma.np,1);
gamx_cheb  = cell(Gamma.np, 1);
gamy_cheb  = cell(Gamma.np, 1);
gam1x_cheb = cell(Gamma.np, 1);
gam1y_cheb = cell(Gamma.np, 1);
for k = 1:Gamma.np
    gamx_cheb{k}  = legvals2chebvals(real(Gamma.x{k}(:,1)));
    gamy_cheb{k}  = legvals2chebvals(real(Gamma.x{k}(:,2)));
    gam1x_cheb{k} = legvals2chebvals(real(Gamma1.x{k}(:,1)));
    gam1y_cheb{k} = legvals2chebvals(real(Gamma1.x{k}(:,2)));
    % Upsample
    gamx_cheb{k}  = bary(xsem, gamx_cheb{k});
    gamy_cheb{k}  = bary(xsem, gamy_cheb{k});
    gam1x_cheb{k} = bary(xsem, gam1x_cheb{k});
    gam1y_cheb{k} = bary(xsem, gam1y_cheb{k});
    xx{k} = t.*gam1x_cheb{k}.' + (1-t).*gamx_cheb{k}.';
    yy{k} = t.*gam1y_cheb{k}.' + (1-t).*gamy_cheb{k}.';
end
dom = cell2struct([xx yy], {'x','y'}, 2);
toc

tic
tfun_un  = cell(length(dom), 1);
tfun  = cell(length(dom), 1);
tfun1 = cell(length(dom), 1);
yfun1d = cell(length(dom), 1);
xfun = cell(length(dom), 1);
yfun = cell(length(dom), 1);
dxfun_dt = cell(length(dom), 1);
dxfun_dr = cell(length(dom), 1);
dyfun_dt = cell(length(dom), 1);
dyfun_dr = cell(length(dom), 1);
for k = 1:length(dom)
    tfun_un{k}  = chebfun(gamx_cheb{k});
    tfun{k}  = chebfun(gamx_cheb{k}, Gamma.breaks(k:k+1));
    tfun1{k} = chebfun(gam1x_cheb{k}, Gamma.breaks(k:k+1));
    xfun{k} = chebfun2(dom(k).x);
    yfun{k} = chebfun2(dom(k).y);
    %xfun{k} = chebfun2(dom(k).x, [Gamma.breaks(k:k+1) -1 1]);
    %yfun{k} = chebfun2(dom(k).y, [Gamma.breaks(k:k+1) -1 1]);
    %scl = 2 / diff(Gamma.breaks(k:k+1));
    scl = 1;
    dxfun_dt{k} = diff(xfun{k}, 1, 2) * scl;
    dxfun_dr{k} = diff(xfun{k}, 1, 1);
    dyfun_dt{k} = diff(yfun{k}, 1, 2) * scl;
    dyfun_dr{k} = diff(yfun{k}, 1, 1);
end
toc

%%
plot(diff(tfun{10}))
hold on
plot(diff(tfun{9}))
shg

%%
fun = xfun;
s = 1;
w = 60;
for k = s:s+w
    surf(xfun{k},yfun{k},fun{k}), hold on
end
plot(Gamma), plot(Gamma1)
view(2)
axis([minmax([dom(s:s+w).x]) minmax([dom(s:s+w).y])])
shg

%% Construct boxes

tic

rGamma1 = refine(Gamma1); % Refined boundary for isinterior() calls.

nb = 16;
% blend_width = 2*nb;
blend_width = 50;
step = util.makeMask(1000, blend_width);

boxes = struct();
%boxes(1).domain   = dom_global;
%boxes(1).domain   = [0 1 0 1];
boxes(1).domain   = [0.25 0.75 0.25 0.75];
boxes(1).level    = 0;
boxes(1).height   = 0;
boxes(1).id       = 1;
boxes(1).parent   = 0;
boxes(1).children = [];
boxes(1).coeffs   = [];
boxes(1).col      = 1;
boxes(1).row      = 1;
elems = {1:Gamma.np};

func = @(x,y) sin(8*pi*x).*sin(8*pi*y);

tol = 1e-11;
m = 2*nb; % Sample at m equispaced points to test error
[xx0, yy0] = chebpts2(nb, nb, [0 1 0 1]);
[xxx0, yyy0] = meshgrid(linspace(0, 1, m));

% Note: the length changes at each iteration here
id = 1;
while ( id <= length(boxes) )

    ncheck = length(boxes)-id+1;
    xx  = zeros(nb, nb, ncheck);
    yy  = zeros(nb, nb, ncheck);
    xxx = zeros(m, m, ncheck);
    yyy = zeros(m, m, ncheck);

    j = 1;
    for k = id:length(boxes)
        % Assemble point list to call isinterior() once
        domk = boxes(k).domain;
        sclx = diff(domk(1:2));
        scly = diff(domk(3:4));
        xx(:,:,j) = sclx*xx0 + domk(1); xxx(:,:,j) = sclx*xxx0 + domk(1);
        yy(:,:,j) = scly*yy0 + domk(3); yyy(:,:,j) = scly*yyy0 + domk(3);
        j = j + 1;
    end
    inGamma_cheb  = isinterior(Gamma,  xx, yy);    inGamma_lin  = isinterior(Gamma,  xxx, yyy);
    inGamma1_cheb = isinterior(rGamma1, xx, yy);   inGamma1_lin = isinterior(rGamma1, xxx, yyy);
    inStrip_cheb  = inGamma_cheb & ~inGamma1_cheb; inStrip_lin  = inGamma_lin & ~inGamma1_lin;
    vals = zeros(nb, nb, ncheck);
    F    = zeros(m, m, ncheck);
    vals(inGamma_cheb) = f(xx(inGamma_cheb), yy(inGamma_cheb));
    F(inGamma_lin)     = f(xxx(inGamma_lin), yyy(inGamma_lin));
    %vals = func(xx,yy);
    %F = func(xxx,yyy);

    % Now evaluate at the strip points
    j = 1;
    for k = id:length(boxes)
        
%         vals_j = zeros(n);
%         vals_j(inGamma_cheb) = f(xx(inGamma_cheb(:,:,j)), yy(inGamma_cheb(:,:,j)));
%         vals

        inStrip_cheb_j = inStrip_cheb(:,:,j);
        inStrip_lin_j  = inStrip_lin(:,:,j);
        xx_j = xx(:,:,j);     xx_sj = xx_j(inStrip_cheb_j);
        yy_j = yy(:,:,j);     yy_sj = yy_j(inStrip_cheb_j);
        vals_j = vals(:,:,j); vals_sj = vals_j(inStrip_cheb_j);
        xxx_j = xxx(:,:,j);   xxx_sj = xxx_j(inStrip_lin_j);
        yyy_j = yyy(:,:,j);   yyy_sj = yyy_j(inStrip_lin_j);
        F_j = F(:,:,j);       F_sj = F_j(inStrip_lin_j);
        rr_sj  = zeros(size(xx_sj));
        rrr_sj = zeros(size(xxx_sj));
        found_cheb = false(size(xx_sj));
        found_lin  = false(size(xxx_sj));

        for e = elems{k}
            % Assume the element list is already filtered
            % Filter out the points that are not in the bounding box of this element
            pad = 0.0001;
            box = [min(dom(e).x(:))-pad max(dom(e).x(:))+pad;
                   min(dom(e).y(:))-pad max(dom(e).y(:))+pad];
            %bnd = 1 + 1e-10;
            bnd = 1 + 1e-3;
            
            %init = [mean(Gamma.breaks(e:e+1)) 0];
            init = [0 0];

            testIdx = find(xx_sj > box(1,1) & xx_sj < box(1,2) & ...
                           yy_sj > box(2,1) & yy_sj < box(2,2));
            if ( ~isempty(testIdx) )
                [tt, rr, converged] = cart2curv_element(xfun{e}, yfun{e}, xx_sj(testIdx), yy_sj(testIdx), ...
                    init, dxfun_dt{e}, dxfun_dr{e}, dyfun_dt{e}, dyfun_dr{e});
                goodIdx = ~found_cheb(testIdx) & converged & abs(tt) <= bnd & abs(rr) <= bnd;
                %goodIdx = ~found_cheb(testIdx) & converged;
                rr_sj(testIdx(goodIdx)) = rr(goodIdx);
                %if ( any(goodIdx) )
                %    vals_sj(testIdx(goodIdx)) = func(xfun{e}(tt(goodIdx),rr(goodIdx)), yfun{e}(tt(goodIdx),rr(goodIdx)));
                %end
                found_cheb(testIdx(goodIdx)) = true;
            end

            testIdx = find(xxx_sj > box(1,1) & xxx_sj < box(1,2) & ...
                           yyy_sj > box(2,1) & yyy_sj < box(2,2));
            if ( ~isempty(testIdx) )
                [ttt, rrr, converged] = cart2curv_element(xfun{e}, yfun{e}, xxx_sj(testIdx), yyy_sj(testIdx), ...
                    init, dxfun_dt{e}, dxfun_dr{e}, dyfun_dt{e}, dyfun_dr{e});
                goodIdx = ~found_lin(testIdx) & converged & abs(ttt) <= bnd & abs(rrr) <= bnd;
                %goodIdx = ~found_lin(testIdx) & converged;
                rrr_sj(testIdx(goodIdx)) = rrr(goodIdx);
                %if ( any(goodIdx) )
                %    F_sj(testIdx(goodIdx)) = func(xfun{e}(ttt(goodIdx),rrr(goodIdx)), yfun{e}(ttt(goodIdx),rrr(goodIdx)));
                %end
                found_lin(testIdx(goodIdx)) = true;
            end
 
        end
        
        %if ( all(abs(boxes(k).domain - [-0.9 -0.825 -0.75 -0.675]) < 1e-11) )
%         if ( all(abs(boxes(k).domain - [-0.9 -0.75 -0.75 -0.6]) < 1e-4) )
%         if ( all(abs(boxes(k).domain - [0.5625 0.567187 -0.76875 -0.76406]) < 1e-4) )
%         if ( all(abs(boxes(k).domain - [0.78125 0.8125 0.46875 0.5]) < 1e-4) )
%             keyboard
%         end
        
%         if ( any(abs(rr_sj) > 1e-10) )
%             keyboard
%         end

        vals_sj = step((rr_sj+1)/2) .* vals_sj;
        F_sj = step((rrr_sj+1)/2) .* F_sj;
        
        if ( ~all(found_cheb) )
            keyboard
        end
        if ( ~all(found_lin) )
            keyboard
        end
        vals_sj(~found_cheb) = 0;
        F_sj(~found_lin) = 0;

        vals_j(inStrip_cheb_j) = vals_sj; % This is like: vals = f(xx,yy);
        vals(:,:,j) = vals_j;
        F_j(inStrip_lin_j) = F_sj;       % This is like: F = f(xxx,yyy)
        F(:,:,j) = F_j;
        coeffs = treefun2.vals2coeffs(vals_j);
        G = treefun2.coeffs2refvals(coeffs);
        err = norm(F_j(:) - G(:), inf);
        resolved = ( err < tol );
        err

        if ( resolved || k >= 1000 )
            boxes(k).coeffs = coeffs;
            boxes(k).height = 0;
        else
            % Split into four child boxes
            domp = boxes(k).domain;
            xmid = mean(domp(1:2));
            ymid = mean(domp(3:4));
            parent = boxes(k);
            
            ex = reshape([dom(elems{k}).x], psem^2, length(elems{k}));
            ey = reshape([dom(elems{k}).y], psem^2, length(elems{k}));
            edom = [min(ex); max(ex); min(ey); max(ey)].';

            child1 = struct();
            child1.domain   = [domp(1) xmid domp(3) ymid];
            child1.id       = length(boxes)+1;
            child1.parent   = k;
            child1.children = [];
            child1.level    = parent.level+1;
            child1.height   = 0;
            child1.coeffs   = [];
            child1.col      = 2*(parent.col-1) + 1;
            child1.row      = 2*(parent.row-1) + 1;
            boxes(end+1) = child1;
            overlapIdx = edom(:,1) < child1.domain(2) & edom(:,2) > child1.domain(1) & ...
                         edom(:,3) < child1.domain(4) & edom(:,4) > child1.domain(3);
            elems{end+1} = elems{k}(overlapIdx);

            child2 = struct();
            child2.domain   = [xmid domp(2) domp(3) ymid];
            child2.id       = length(boxes)+1;
            child2.parent   = k;
            child2.children = [];
            child2.level    = parent.level+1;
            child2.height   = 0;
            child2.coeffs   = [];
            child2.col      = 2*(parent.col-1) + 2;
            child2.row      = 2*(parent.row-1) + 1;
            boxes(end+1) = child2;
            overlapIdx = edom(:,1) < child2.domain(2) & edom(:,2) > child2.domain(1) & ...
                         edom(:,3) < child2.domain(4) & edom(:,4) > child2.domain(3);
            elems{end+1} = elems{k}(overlapIdx);

            child3 = struct();
            child3.domain   = [domp(1) xmid ymid domp(4)];
            child3.id       = length(boxes)+1;
            child3.parent   = k;
            child3.children = [];
            child3.level    = parent.level+1;
            child3.height   = 0;
            child3.coeffs   = [];
            child3.col      = 2*(parent.col-1) + 1;
            child3.row      = 2*(parent.row-1) + 2;
            boxes(end+1) = child3;
            overlapIdx = edom(:,1) < child3.domain(2) & edom(:,2) > child3.domain(1) & ...
                         edom(:,3) < child3.domain(4) & edom(:,4) > child3.domain(3);
            elems{end+1} = elems{k}(overlapIdx);

            child4 = struct();
            child4.domain   = [xmid domp(2) ymid domp(4)];
            child4.id       = length(boxes)+1;
            child4.parent   = k;
            child4.children = [];
            child4.level    = parent.level+1;
            child4.height   = 0;
            child4.coeffs   = [];
            child4.col      = 2*(parent.col-1) + 2;
            child4.row      = 2*(parent.row-1) + 2;
            boxes(end+1) = child4;
            overlapIdx = edom(:,1) < child4.domain(2) & edom(:,2) > child4.domain(1) & ...
                         edom(:,3) < child4.domain(4) & edom(:,4) > child4.domain(3);
            elems{end+1} = elems{k}(overlapIdx);

            boxes(k).children = [child1.id, child2.id, ...
                                 child3.id, child4.id];

            %f.boxes(id).children = [child3.id, child4.id, ...
            %                        child2.id, child1.id];
            
            boxes(k).coeffs = [];            
            boxes(k).height = 1;
        end

        j = j + 1;
    end
    
    id = id + ncheck;    
end

toc

% Convert to treefun
tf = treefun2(boxes);

%%
e = 9;
ebdy_x = linspace(real(Gamma1.zbreaks(e)), real(Gamma.zbreaks(e)), 20).';
ebdy_y = linspace(imag(Gamma1.zbreaks(e)), imag(Gamma.zbreaks(e)), 20).';
%init = [-ones(20,1) linspace(1, -1, 20).'];
init = zeros(101,2);
tt1 = linspace(-1,1,101).';
rr1 = zeros(101,1);
xp = xfun{e}(tt1,rr1);
yp = yfun{e}(tt1,rr1);
[tt, rr, flag] = cart2curv_element(xfun{e}, yfun{e}, xp, yp, ...
    init, dxfun_dt{e}, dxfun_dr{e}, dyfun_dt{e}, dyfun_dr{e});

e = e+1;
xp = xfun{e}(tt1,rr1);
yp = yfun{e}(tt1,rr1);
[tt, rr, flag] = cart2curv_element(xfun{e}, yfun{e}, xp, yp, ...
    init, dxfun_dt{e}, dxfun_dr{e}, dyfun_dt{e}, dyfun_dr{e});

% [tt, rr, flag] = cart2curv_element2(xfun{e}, yfun{e}, 0.62, 0.34, ...
%     [0 0], dxfun_dt{e}, dxfun_dr{e}, dyfun_dt{e}, dyfun_dr{e});
%[tt, rr, converged] = cart2curv_element(xfun{e}, yfun{e}, ebdy_x, ebdy_y, ...
%    init, dxfun_dt{e}, dxfun_dr{e}, dyfun_dt{e}, dyfun_dr{e});
%[tt1, rr1, converged1] = cart2curv_element(xfun{e-1}, yfun{e-1}, ebdy_x, ebdy_y, ...
%    [-1 1], dxfun_dt{e-1}, dxfun_dr{e-1}, dyfun_dt{e-1}, dyfun_dr{e-1});

