close all
fprintf('\n');

%% Define Gamma
n = 16;          % Number of nodes per panel
m = n;           % Number of radial Chebyshev nodes in strip
n_re = 2*n;      % Resample original curve to use n_re per panel
nsem = 2*n;      % Number of nodes for SEM
n_box = n;       % Each box uses (n_box x n_box) Chebyshev nodes
m_box = 2*n_box; % Box error is computed on (m_box x m_box) equispaced points

Gamma = Boundary.circle(n, 'quadrature', 'panel');
%Gamma = Boundary.ellipse(n, 5, 1, 'quadrature', 'panel');
%Gamma = Boundary.squircle(n, 10, 1, 1, 'quadrature', 'panel');
%Gamma = refine(Gamma, 1);
%Gamma = Boundary.star(n, 'quadrature', 'panel');
%Gamma = refine(Gamma);
%Gamma = Boundary.multiscale_circle(n, 'quadrature', 'panel');
%Gamma = Boundary.wavepacket(n, 'quadrature', 'panel');
%Gamma = Boundary.c_shape(n, 'quadrature', 'panel');
%Gamma = refine(Gamma);

tic
Gamma_re = resample(Gamma, n_re);
fprintf('   Resampling Gamma ................... %.6fs\n', toc);

gamxy = cell2mat(Gamma_re.x);
gamx  = gamxy(:,1);
gamy  = gamxy(:,2);

scale = 1.2;
dom_global = boundingbox(Gamma, scale);
dom_global([1 3]) = min(dom_global);
dom_global([2 4]) = max(dom_global);

syms x y
gaussian = @(x,y,cx,cy,K) 1/K*exp(-K^2*((x-cx).^2 + (y-cy).^2));
%sol = @(x,y) 0.5*x.^3 + 0.5*y.^3;% + ...
%sol = @(x,y) -cos(30*x).*exp(sin(x)).*sin(30*y);% + ...
sol = @(x,y) -cos(5*x).*sin(5*y);% + ...
             %gaussian(x,y,0,0,10) + ...
             %gaussian(x,y,0,1-1e-1,100);% + ...
             %gaussian(x,y,0,1-1e-2,1000) + ...
             %gaussian(x,y,0,1-1e-3,10000);
f = matlabFunction(simplify(laplacian(sol(x,y))));
g = @(x,y) sol(x,y);

%% Define the fake interior boundary, Gamma'
tic
beta = 4;
bleed = 10;
width = 0.5;
[x, y] = smoothStrip2(Gamma_re, n_re, beta, bleed, width);
xleg = chebvals2legvals(reshape(x, n_re, Gamma.np)); xleg = xleg(:);
yleg = chebvals2legvals(reshape(y, n_re, Gamma.np)); yleg = yleg(:);
z1 = mat2cell(xleg + 1i*yleg, repmat(n_re, Gamma.np, 1), 1);
Gamma1 = Boundary(z1);

% z = cell(Gamma_re.np, 1);
% width = 0.2;
% for k = 1:Gamma_re.np
%     nx = Gamma_re.normal{k};
%     nx = nx(:,1)+nx(:,2)*1i;
%     nx = nx./abs(nx);
%     z{k} = Gamma_re.z{k} - width*nx;
% end
% Gamma1 = Boundary(z);

% Build the strip grid
xsem = chebpts(n_re);
t = chebpts(n_re, [0 1]);
xx = cell(Gamma_re.np,1);
yy = cell(Gamma_re.np,1);
gamx_cheb  = cell(Gamma_re.np, 1);
gamy_cheb  = cell(Gamma_re.np, 1);
gam1x_cheb = cell(Gamma_re.np, 1);
gam1y_cheb = cell(Gamma_re.np, 1);
gam_nx_cheb = cell(Gamma_re.np, 1);
gam1_nx_cheb = cell(Gamma_re.np, 1);
LV2CV = legvals2chebvals(eye(n_re));
for k = 1:Gamma.np
    gamx_cheb{k}  = LV2CV * real(Gamma_re.x{k}(:,1));
    gamy_cheb{k}  = LV2CV * real(Gamma_re.x{k}(:,2));
    gam1x_cheb{k} = LV2CV * real(Gamma1.x{k}(:,1));
    gam1y_cheb{k} = LV2CV * real(Gamma1.x{k}(:,2));
    gam_nx_cheb{k} = LV2CV * Gamma_re.normal{k};
    gam1_nx_cheb{k} = LV2CV * Gamma1.normal{k};
    xx{k} = t.*gam1x_cheb{k}.' + (1-t).*gamx_cheb{k}.';
    yy{k} = t.*gam1y_cheb{k}.' + (1-t).*gamy_cheb{k}.';
end
dom = cell2struct([xx yy], {'x','y'}, 2);

fprintf('   Constructing Gamma'' ................ %.6fs\n', toc);

%% Construct boxes
tic

rGamma1 = refine(Gamma1); % Refined boundary for isinterior() calls.
%rGamma1 = Gamma1; % Refined boundary for isinterior() calls.

% blend_width = 2*n_box;
blend_width = 30;
step = util.blender(blend_width, 'pswf', [0 1]);
step_cfs = step.coeffs;
step = @(x) chebtech.clenshaw(x, step_cfs);

% Allocate space for many boxes
maxBoxes = 20000;
domains  = zeros(4, maxBoxes);
levels   = zeros(maxBoxes, 1);
heights  = zeros(maxBoxes, 1);
ids      = zeros(maxBoxes, 1);
parents  = zeros(maxBoxes, 1);
children = zeros(4, maxBoxes);
coeffs   = cell(maxBoxes, 1);
cols     = zeros(maxBoxes, 1);
rows     = zeros(maxBoxes, 1);
elems    = cell(maxBoxes, 1);
boxInStrip  = true(maxBoxes, 1);

% Assign values for the root box
nboxes = 1;
domains(:,1) = dom_global;
levels(1)    = 0;
heights(1)   = 0;
ids(1)       = 1;
parents(1)   = 0;
cols(1)      = 1;
rows(1)      = 1;
elems{1} = 1:Gamma.np;
boxInStrip(1) = true;

nthreads = 8;
tol = 1e-6;
[xx0, yy0] = chebpts2(n_box, n_box, [0 1 0 1]);
[xxx0, yyy0] = meshgrid(linspace(0, 1, m_box));
[xcheb, ~, vcheb] = chebpts(n_re);
getElemCoordinates = util.getElemCoordinates_(n_re);

% Note: the length changes at each iteration here
id = 1;
while ( id <= nboxes )
    id

    ncheck = nboxes-id+1;
    xx  = zeros(n_box, n_box, ncheck);
    yy  = zeros(n_box, n_box, ncheck);
    xxx = zeros(m_box, m_box, ncheck);
    yyy = zeros(m_box, m_box, ncheck);

    j = 1;
    for k = id:nboxes
        % Assemble point list to call isinterior() once
        domk = domains(:,k);
        sclx = diff(domk(1:2));
        scly = diff(domk(3:4));
        xx(:,:,j) = sclx*xx0 + domk(1); xxx(:,:,j) = sclx*xxx0 + domk(1);
        yy(:,:,j) = scly*yy0 + domk(3); yyy(:,:,j) = scly*yyy0 + domk(3);
        j = j + 1;
    end
    inGamma_cheb  = isinterior(Gamma,  xx, yy);    inGamma_lin  = isinterior(Gamma,  xxx, yyy);
    inGamma1_cheb = isinterior(rGamma1, xx, yy);   inGamma1_lin = isinterior(rGamma1, xxx, yyy);
    inStrip_cheb  = inGamma_cheb & ~inGamma1_cheb; inStrip_lin  = inGamma_lin & ~inGamma1_lin;
    vals = zeros(n_box, n_box, ncheck);
    F    = zeros(m_box, m_box, ncheck);
    vals(inGamma_cheb) = f(xx(inGamma_cheb), yy(inGamma_cheb));
    F(inGamma_lin)     = f(xxx(inGamma_lin), yyy(inGamma_lin));

    % Now evaluate at the strip points
    j = 1;
    nboxes_old = nboxes;
    for k = id:nboxes_old

        inStrip_cheb_j = inStrip_cheb(:,:,j);
        inStrip_lin_j  = inStrip_lin(:,:,j);
        boxInStrip(k) = any(inStrip_lin_j, 'all');
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
            pad = 10*tol;
            box = [min(dom(e).x(:))-1e-12 max(dom(e).x(:))+1e-12;
                   min(dom(e).y(:))-1e-12 max(dom(e).y(:))+1e-12];
            aL = [dom(e).x(1,1)   dom(e).y(1,1)];
            bL = [dom(e).x(end,1) dom(e).y(end,1)];
            aR = [dom(e).x(1,end)   dom(e).y(1,end)];
            bR = [dom(e).x(end,end) dom(e).y(end,end)];
            
            nx_j = gam_nx_cheb{e}(1,:);
            perpx = -nx_j(2);
            perpy =  nx_j(1);
            
            inBoundingBox = xx_sj > box(1,1) & xx_sj < box(1,2) & ...
                            yy_sj > box(2,1) & yy_sj < box(2,2);
            inElem = inBoundingBox & ~found_cheb;
            if ( any(inElem) )
                nx = gam_nx_cheb{e}(:,1);
                ny = gam_nx_cheb{e}(:,2);
                [~, rr, found] = getElemCoordinates(xx_sj(inElem), ...
                    yy_sj(inElem), gamx_cheb{e}, gamy_cheb{e}, ...
                    gam1x_cheb{e}, gam1y_cheb{e}, nx, ny, xcheb, vcheb, nthreads);
                if ( any(isnan(rr(found))) )
                    keyboard
                end
                found(abs(rr) > 1) = false;
                inElem(inElem) = found;
                rr_sj(inElem) = rr(found);
                found_cheb(inElem) = true;
            end

            inBoundingBox = xxx_sj > box(1,1) & xxx_sj < box(1,2) & ...
                            yyy_sj > box(2,1) & yyy_sj < box(2,2);
            inElem = inBoundingBox & ~found_lin;
            if ( any(inElem) )
                nx = gam_nx_cheb{e}(:,1);
                ny = gam_nx_cheb{e}(:,2);
                [~, rrr, found] = getElemCoordinates(xxx_sj(inElem), ...
                    yyy_sj(inElem), gamx_cheb{e}, gamy_cheb{e}, ...
                    gam1x_cheb{e}, gam1y_cheb{e}, nx, ny, xcheb, vcheb, nthreads);
                found(abs(rrr) > 1) = false;
                inElem(inElem) = found;
                rrr_sj(inElem) = rrr(found);
                found_lin(inElem) = true;
            end
 
        end

        ss = step([rr_sj ; rrr_sj]);
        vals_sj = ss(1:length(rr_sj)) .* vals_sj;
        F_sj = ss(length(rr_sj)+1:end) .* F_sj;
        
        if ( any(abs(vals_sj) > 1e10) )
            keyboard
        end

        vals_sj(~found_cheb) = 0;
        F_sj(~found_lin) = 0;

        vals_j(inStrip_cheb_j) = vals_sj; % This is like: vals = f(xx,yy);
        vals(:,:,j) = vals_j;
        F_j(inStrip_lin_j) = F_sj;        % This is like: F = f(xxx,yyy)
        F(:,:,j) = F_j;
        %coeffs_j = treefun2.vals2coeffs(vals_j);
        %G = treefun2.coeffs2refvals(coeffs_j);
        G = treefun2.chebvals2refvals(vals_j);
        err = norm(F_j(:) - G(:), inf);
        
        %coeffs_j = treefun2.vals2coeffs(vals_j);
        %Ex = sum(abs(coeffs_j(end-1:end,:)), 'all') / (2*n);
        %Ey = sum(abs(coeffs_j(:,end-1:end)), 'all') / (2*n);
        %err_cfs = (Ex + Ey) / 2;
        
        %vscale = max( abs( vals_j(:) ) );
        %resolved = ( err_cfs < tol );
        resolved = ( err < tol );
        
        if ( isnan(err) )
            keyboard
        end

        if ( resolved )
            coeffs{k} = treefun2.vals2coeffs(vals_j);
            heights(k) = 0;
        else
            % Split into four child boxes
            domp = domains(:,k);
            xmid = mean(domp(1:2));
            ymid = mean(domp(3:4));
            
            ex = reshape([dom(elems{k}).x], n_re^2, length(elems{k}));
            ey = reshape([dom(elems{k}).y], n_re^2, length(elems{k}));
            edom = [min(ex); max(ex); min(ey); max(ey)].';

            cid1 = nboxes+1;
            domains(:,cid1) = [domp(1) xmid domp(3) ymid];
            ids(cid1)       = cid1;
            parents(cid1)   = k;
            levels(cid1)    = levels(k)+1;
            heights(cid1)   = 0;
            cols(cid1)      = 2*(cols(k)-1) + 1;
            rows(cid1)      = 2*(rows(k)-1) + 1;
            overlapIdx = edom(:,1) < domains(2,cid1) & edom(:,2) > domains(1,cid1) & ...
                         edom(:,3) < domains(4,cid1) & edom(:,4) > domains(3,cid1);
            elems{cid1} = elems{k}(overlapIdx);
            if ( isempty(elems{cid1}) )
                elems{cid1} = [];
            end

            cid2 = nboxes+2;
            domains(:,cid2) = [xmid domp(2) domp(3) ymid];
            ids(cid2)       = cid2;
            parents(cid2)   = k;
            levels(cid2)    = levels(k)+1;
            heights(cid2)   = 0;
            cols(cid2)      = 2*(cols(k)-1) + 2;
            rows(cid2)      = 2*(rows(k)-1) + 1;
            overlapIdx = edom(:,1) < domains(2,cid2) & edom(:,2) > domains(1,cid2) & ...
                         edom(:,3) < domains(4,cid2) & edom(:,4) > domains(3,cid2);
            elems{cid2} = elems{k}(overlapIdx);
            if ( isempty(elems{cid2}) )
                elems{cid2} = [];
            end

            cid3 = nboxes+3;
            domains(:,cid3) = [domp(1) xmid ymid domp(4)];
            ids(cid3)       = cid3;
            parents(cid3)   = k;
            levels(cid3)    = levels(k)+1;
            heights(cid3)   = 0;
            cols(cid3)      = 2*(cols(k)-1) + 1;
            rows(cid3)      = 2*(rows(k)-1) + 2;
            overlapIdx = edom(:,1) < domains(2,cid3) & edom(:,2) > domains(1,cid3) & ...
                         edom(:,3) < domains(4,cid3) & edom(:,4) > domains(3,cid3);
            elems{cid3} = elems{k}(overlapIdx);
            if ( isempty(elems{cid3}) )
                elems{cid3} = [];
            end

            cid4 = nboxes+4;
            domains(:,cid4) = [xmid domp(2) ymid domp(4)];
            ids(cid4)       = cid4;
            parents(cid4)   = k;
            levels(cid4)    = levels(k)+1;
            heights(cid4)   = 0;
            cols(cid4)      = 2*(cols(k)-1) + 2;
            rows(cid4)      = 2*(rows(k)-1) + 2;
            overlapIdx = edom(:,1) < domains(2,cid4) & edom(:,2) > domains(1,cid4) & ...
                         edom(:,3) < domains(4,cid4) & edom(:,4) > domains(3,cid4);
            elems{cid4} = elems{k}(overlapIdx);
            if ( isempty(elems{cid4}) )
                elems{cid4} = [];
            end

            if ( ~boxInStrip(k) )
                boxInStrip([cid1 cid2 cid3 cid4]) = false;
            end
            children(:,k) = [cid1 cid2 cid3 cid4];
            heights(k) = 1;
            nboxes = nboxes + 4;
        end

        j = j + 1;
    end
    
    id = id + ncheck;
    
    if ( nboxes > maxBoxes )
        keyboard
    end
end

% Do a cumulative sum in reverse to correct the heights
for k = nboxes:-1:1
    if ( heights(k) ~= 0 )
        heights(k) = 1 + max(heights(children(:,k)));
    end
end

domains  = domains(:,1:nboxes);
levels   = levels(1:nboxes);
heights  = heights(1:nboxes);
ids      = ids(1:nboxes);
parents  = parents(1:nboxes);
children = children(:,1:nboxes);
coeffs   = coeffs(1:nboxes);
cols     = cols(1:nboxes);
rows     = rows(1:nboxes);

% Convert to treefun (this will balance the tree)
tf = treefun2_arrays(domains, levels, heights, ids, parents, children, ...
   coeffs, cols, rows);

fprintf('   Constructing quadtree .............. %.6fs\n', toc);

%% (1) Solve the bulk problem
%
%    lap(u_bulk) = phi*f     in B
tic

u_bulk = treefun2.poisson(-tf);
%u_bulk_true = chebfun2(sol, dom_global);
%u_bulk_true = treefun2(@(x,y) sol(x,y), u_bulk);
%u_bulk = u_bulk_true;

fprintf('   Box code ........................... %.6fs\n', toc);

%% (2) Solve the strip problem
%
%    lap(u_strip) = f  in S
%               u = 0  on Gamma & Gamma'
tic

% Build the strip grid
xsem = chebpts(nsem);
t = chebpts(nsem, [0 1]);
xx = cell(Gamma_re.np,1);
yy = cell(Gamma_re.np,1);
gamx_cheb_resamp  = cell(Gamma_re.np, 1);
gamy_cheb_resamp  = cell(Gamma_re.np, 1);
gam1x_cheb_resamp = cell(Gamma_re.np, 1);
gam1y_cheb_resamp = cell(Gamma_re.np, 1);
gam_nx_cheb_resamp = cell(Gamma_re.np, 1);
gam1_nx_cheb_resamp = cell(Gamma_re.np, 1);
for k = 1:Gamma.np
    gamx_cheb_resamp{k}  = bary(xsem, gamx_cheb{k});
    gamy_cheb_resamp{k}  = bary(xsem, gamy_cheb{k});
    gam1x_cheb_resamp{k} = bary(xsem, gam1x_cheb{k});
    gam1y_cheb_resamp{k} = bary(xsem, gam1y_cheb{k});
    gam_nx_cheb_resamp{k} = bary(xsem, gam_nx_cheb{k});
    gam1_nx_cheb_resamp{k} = bary(xsem, gam1_nx_cheb{k});
    xx{k} = t.*gam1x_cheb_resamp{k}.' + (1-t).*gamx_cheb_resamp{k}.';
    yy{k} = t.*gam1y_cheb_resamp{k}.' + (1-t).*gamy_cheb_resamp{k}.';
end
dom_sem = cell2struct([xx yy], {'x','y'}, 2);

% dom_sem = dom;
% nsem = n;

S = StripSolver(dom_sem, f);
build(S);

outerIdx = false(size(S.patches{1}.xy, 1), 1);
outerIdx((1:nsem-2).' + (nsem-2)*(0:2:2*length(dom_sem)-1)) = true;
innerIdx = ~outerIdx;
innerx = S.patches{1}.xy(innerIdx, 1);
innery = S.patches{1}.xy(innerIdx, 2);
outerx = S.patches{1}.xy(outerIdx, 1);
outery = S.patches{1}.xy(outerIdx, 2);

% This will call treefun2/feval
bc = u_bulk(S.patches{1}.xy(:,1), S.patches{1}.xy(:,2));
%bc = u_bulk_true(S.patches{1}.xy(:,1), S.patches{1}.xy(:,2));
%bc_true = u_bulk_true(S.patches{1}.xy(:,1), S.patches{1}.xy(:,2));
%bc(outerIdx) = bc_true(outerIdx);
%bc(outerIdx) = bc(innerIdx);
%bc(outerIdx) = g(outerx, outery);

ubc = reshape(bc(outerIdx), nsem-2, length(dom));
%ubc = util.modchebvals2legvals(ubc);
ubc = util.modchebvals2chebvals(ubc);
x_re = legpts(n_re);
ubc = bary(x_re, ubc);
ubc = ubc(:);

u_strip = S \ bc;

% u_strip_true = cell(length(dom_sem), 1);
% for k = 1:length(dom_sem)
%     u_strip_true{k} = sol(dom_sem(k).x, dom_sem(k).y);
% end
% u_strip = u_strip_true;
% gamxy_re = cell2mat(Gamma_re.x);
% gamx_re  = gamxy_re(:,1);
% gamy_re  = gamxy_re(:,2);
% ubc = sol(gamx_re,gamy_re);

fprintf('   Strip solve ........................ %.6fs\n', toc);

%% (3) Solve the Neumann glue problem
%
%     lap(u_glue) = 0             in B \ Gamma'
%        [u_glue] = 0             on Gamma'
%    [du_glue/dn] = -du_strip/dn  on Gamma'
tic

% Compute Neumann jump on Gamma'
normal1 = cell2mat(Gamma1.normal);
nx1 = normal1(:,1);
ny1 = normal1(:,2);
gam1x = cell2mat(Gamma1.x); gam1x = gam1x(:,1);
gam1y = cell2mat(Gamma1.x); gam1y = gam1y(:,2);

u_bulk_dn = feval(diff(u_bulk,1,2), gam1x, gam1y) .* nx1 + ...
            feval(diff(u_bulk,1,1), gam1x, gam1y) .* ny1;

u_strip_dn_all = S.patches{1}.D2N * [bc ; 1];
u_strip_dn = reshape(u_strip_dn_all(innerIdx), nsem-2, length(dom));
%u_strip_dn = util.modchebvals2legvals(u_strip_dn);
u_strip_dn = util.modchebvals2chebvals(u_strip_dn);
x_re = legpts(n_re);
u_strip_dn = bary(x_re, u_strip_dn);
u_strip_dn = u_strip_dn(:);

% syms x y
% u_strip_dx = matlabFunction(diff(sol(x,y), 'x'));
% u_strip_dy = matlabFunction(diff(sol(x,y), 'y'));
% u_strip_dn = -feval(u_strip_dx, gam1x, gam1y) .* nx1 + ...
%              -feval(u_strip_dy, gam1x, gam1y) .* ny1;

neu_jump = -(u_strip_dn + u_bulk_dn);
u_glue_i = @(x,y) kernels.laplace.slp(Gamma1, 'density', neu_jump, 'target', [x y], 'closeeval', true, 'side', 'i');
u_glue_e = @(x,y) kernels.laplace.slp(Gamma1, 'density', neu_jump, 'target', [x y], 'closeeval', true, 'side', 'e');

fprintf('   Correcting normal derivative ....... %.6fs\n', toc);

%% (4) Solve the boundary correction problem
%
%    lap(u_bc) = 0                                  in B \ Gamma
%         u_bc = g - u_bulk|_Gamma - u_glue|_Gamma  on Gamma
tic
% 
% gamxy = cell2mat(Gamma.x);
% gamx  = gamxy(:,1);
% gamy  = gamxy(:,2);
% bc = u_bulk(S.patches{1}.xy(:,1), S.patches{1}.xy(:,2));
% ubc = reshape(bc(outerIdx), nsem-2, length(dom));
% ubc = util.modchebvals2chebvals(ubc);
% x_re = legpts(n);
% ubc = bary(x_re, ubc);
% ubc = ubc(:);
% bc = g(gamx,gamy) - ubc - u_glue_e(gamx,gamy);
% chnkr1 = makeChunker(Gamma);
% fkern = @(s,t) chnk.lap2d.kern(s,t,'D');
% D = chunkermat(chnkr1, fkern);
% I = eye(n * Gamma.np);
% sigma = (D - I/2) \ bc;
% %sigma = gmres(D - I/2, bc, [], 1e-14, 50);

gamxy_re = cell2mat(Gamma_re.x);
gamx_re  = gamxy_re(:,1);
gamy_re  = gamxy_re(:,2);
bc = g(gamx_re,gamy_re) - ubc - u_glue_e(gamx_re,gamy_re);
K = kernels.laplace.dlp(Gamma_re);
I = eye(n_re * Gamma_re.np);
sigma = (K - I/2) \ bc;
% sigma = gmres(K - I/2, bc, [], 1e-14, 50);

u_bc = @(x,y) kernels.laplace.dlp(Gamma_re, 'density', sigma, 'target', [x y], 'closeeval', true, 'side', 'i');
%u_bc = @(x,y) chunkerkerneval(chnkr1, fkern, sigma, [x(:).'; y(:).']);

fprintf('   Correcting boundary data ........... %.6fs\n', toc);

%% Add up the solutions

tic
nbulk = 400;
[xx, yy] = meshgrid(linspace(dom_global(1), dom_global(2), nbulk).', ...
                    linspace(dom_global(3), dom_global(4), nbulk).');
%xx = xx(:);
%yy = yy(:);
%[xx, yy] = leafpts(tf);
%xx = [xx{:}]; xx = xx(:);
%yy = [yy{:}]; yy = yy(:);
inGamma  = isinterior(Gamma_re,  xx, yy);
inGamma1 = isinterior(Gamma1, xx, yy);
inStrip  = inGamma & ~inGamma1;
xx_g = xx(inGamma); xx_g1 = xx(inGamma1); xx_s = xx(inStrip);
yy_g = yy(inGamma); yy_g1 = yy(inGamma1); yy_s = yy(inStrip);
%[tt_s, rr_s, globalIdx] = util.assignStripPoints(dom, xx_s, yy_s);
[tt_s, rr_s, globalIdx] = util.assignStripPoints(dom, xx_s, yy_s, ...
    gamx_cheb, gamy_cheb, gam1x_cheb, gam1y_cheb, gam_nx_cheb);
fprintf('   Getting plot points ................ %.6fs\n', toc);

tic
uu = nan(nbulk);
%uu = nan(numel(xx),1);
ub = u_bulk(xx,yy);
uu(inGamma1) = ub(inGamma1) + u_glue_i(xx_g1,yy_g1) + u_bc(xx_g1,yy_g1);
stripvals = zeros(size(xx_s));
for k = 1:length(dom)
    stripvals(globalIdx{k}) = util.bary2d(u_strip{k}, tt_s{k}, rr_s{k});
end
uu(inStrip) = stripvals + u_glue_e(xx_s, yy_s) + u_bc(xx_s, yy_s);
fprintf('   Evaluating solution ................ %.6fs\n', toc);

uu_sol = nan(nbulk);
%uu_sol = nan(numel(xx),1);
uu_sol(inGamma) = sol(xx_g,yy_g);
% uu_sol1 = nan(nbulk);
% uu_sol1(inGamma1) = sol(xx_g1,yy_g1);
% uu_solstrip = nan(nbulk);
% uu_solstrip(inStrip) = sol(xx_s,yy_s);
% 
% ug = nan(size(xx));
% ug(inGamma1) = uu(inGamma1);

max(max(abs(uu - uu_sol)))

uu = real(uu);

figure(1)
%surf(xx, yy, ug - uu_sol1)
%surf(xx, yy, uu_sol1)
%surf(xx, yy, uu)
%surf(xx, yy, uu-uu_sol)
surf(xx, yy, log10(abs(uu-uu_sol)))
shading interp
shg
