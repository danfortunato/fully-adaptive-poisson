function tf = constructTreeIntension(S, f)

tic

Gamma = S.Gamma;
Gamma1 = S.Gamma1;
n_re  = S.n_re;
n_box = S.n_box;
m_box = S.m_box;
dom = S.strip_dom;
gamx_cheb   = S.gamx_cheb;
gamy_cheb   = S.gamy_cheb;
gam1x_cheb  = S.gam1x_cheb;
gam1y_cheb  = S.gam1y_cheb;
gam_nx_cheb = S.gam_nx_cheb;

% Refined boundary for isinterior() calls.
rGamma1 = refine(Gamma1);
%rGamma1 = Gamma1;

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
domains(:,1) = S.domain;
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
% tf = treefun2_arrays(domains, levels, heights, ids, parents, children, ...
%    coeffs, cols, rows);

boxes = struct();
for k = 1:nboxes
    boxes(k).domain = domains(:,k).';
    boxes(k).level = levels(k);
    boxes(k).height = heights(k);
    boxes(k).id = ids(k);
    boxes(k).parent = parents(k);
    boxes(k).children = children(:,k).';
    boxes(k).coeffs = coeffs{k};
    boxes(k).col = cols(k);
    boxes(k).row = rows(k);
end
tf = treefun2(boxes);

fprintf('   Constructing quadtree .............. %.6fs\n', toc);

end