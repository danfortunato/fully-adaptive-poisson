function [tf, tf1] = constructTreeIntension2(S, f)

tic

nthreads = 8;
tol = 1e-8;

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
% rGamma1 = refine(Gamma1);
rGamma1 = Gamma1;

% blend_width = 2*n_box;
blend_width = 30;
step = util.blender(blend_width, 'pswf', [0 1]);
step_cfs = step.coeffs;
step = @(x) chebtech.clenshaw(x, step_cfs);

% Allocate space for many boxes
maxBoxes = 100000;
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

% Flags to tell where the boxes lie in relation to each curve:
%    1 = fully inside
%   -1 = fully outside
%    0 = neither
boxInGamma  = zeros(maxBoxes, 1);
boxInGamma1 = zeros(maxBoxes, 1);

% Make KD-trees for Gamma and Gamma'
kd  = KDTree(cell2mat(Gamma.x));
kd1 = KDTree(cell2mat(Gamma1.x));

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

[xx0, yy0] = chebpts2(n_box, n_box, [0 1 0 1]);
[xxx0, yyy0] = meshgrid(linspace(0, 1, m_box));
[xcheb, ~, vcheb] = chebpts(n_re);
getElemCoordinates = util.getElemCoordinates_(n_re);

% Error model to use for resolve checking: 'vals' or 'coeffs'
errorType = 'coeffs';

% Note: the length changes at each iteration here
id = 1;
while ( id <= nboxes )
    id

    % We will check all the leaves on this level for refinement.
    ncheck = nboxes-id+1;

    % Try to avoid calling isinterior on the dense set of leaf points.
    % We will proceed in three phases.

    %%% Phase 1: If our parent is fully inside or outside of both curves
    %   then we know where we are.
    if ( id > 1 )
        boxInGamma(id:nboxes)  = boxInGamma(parents(id:nboxes));
        boxInGamma1(id:nboxes) = boxInGamma1(parents(id:nboxes));
    end
    boxesToCheck  = ( boxInGamma(id:nboxes)  == 0 );
    boxesToCheck1 = ( boxInGamma1(id:nboxes) == 0 );
    
    %%% Phase 2: If we weren't able to determine where we are from Phase 1,
    %   then we check where the four corners of our box lie. If they are
    %   all on the same side of the curve and no boundary lie near it, then
    %   we know where we are. We use a KD-tree to perform this range query.
    
    % Collect all the corners of the questionable boxes into an array so we
    % can call isinterior() once, on the questionable ones.
    corners_x = domains([1 2 2 1], id:nboxes);
    corners_y = domains([3 3 4 4], id:nboxes);
    cornersInGamma  = isinterior(Gamma,  corners_x(:,boxesToCheck),  corners_y(:,boxesToCheck));
    cornersInGamma1 = isinterior(Gamma1, corners_x(:,boxesToCheck1), corners_y(:,boxesToCheck1));

    j = 1; % Too many indices...
    l = 1;
    l1 = 1;
    for k = id:nboxes
        % We will look for any points in Gamma or Gamma' that lie "near"
        % this box, where "near" means inside of a slightly larger box
        % enclosing this one.
        if ( boxesToCheck(j) || boxesToCheck1(j) )
            % Pad by a multiple of the diagonal length of the box
            pad = 0.1*sqrt(2)*diff(domains(1:2, k));
            box = [domains(1,k)-pad domains(2,k)+pad;
                   domains(3,k)-pad domains(4,k)+pad];
        end
        if ( boxesToCheck(j) )
            % This box's relation to Gamma is still unknown. Check if all
            % the corners lie on the same side and use the KD tree to test
            % for nearness.
            nearGamma = ~isempty(kd.range(box));
            if ( all(cornersInGamma(:,l)) && ~nearGamma )
                boxInGamma(k) = 1;
            elseif ( all(~cornersInGamma(:,l)) && ~nearGamma )
                boxInGamma(k) = -1;
            end
            l = l + 1;
        end
        if ( boxesToCheck1(j) )
            % This box's relation to Gamma is still unknown. Check if all
            % the corners lie on the same side and use the KD tree to test
            % for nearness.
            nearGamma1 = ~isempty(kd1.range(box));
            if ( all(cornersInGamma1(:,l1)) && ~nearGamma1 )
                boxInGamma1(k) = 1;
            elseif ( all(~cornersInGamma1(:,l1)) && ~nearGamma1 )
                boxInGamma1(k) = -1;
            end
            l1 = l1 + 1;
        end
        j = j + 1;
    end
    
    % Update the boxes that still checking.
    boxesToCheck  = ( boxInGamma(id:nboxes)  == 0 );
    boxesToCheck1 = ( boxInGamma1(id:nboxes) == 0 );
    
    %%% Phase 3: We still don't know where this box lies in relation to
    % Gamma and/or Gamma'. We have no choice but to call isinterior() on
    % the dense set of leaf nodes.

    xx  = zeros(n_box, n_box, ncheck);
    yy  = zeros(n_box, n_box, ncheck);
    xxx = zeros(m_box, m_box, ncheck);
    yyy = zeros(m_box, m_box, ncheck);
    inGamma_cheb  = false(n_box, n_box, ncheck);
    inGamma_lin   = false(m_box, m_box, ncheck);
    inGamma1_cheb = false(n_box, n_box, ncheck);
    inGamma1_lin  = false(m_box, m_box, ncheck);

    j = 1;
    for k = id:nboxes
        % Fill in the point locations for the boxes we know.
        % Gamma locations:
        if ( boxInGamma(k) == 1 )
            inGamma_cheb(:,:,j) = true;  inGamma_lin(:,:,j) = true;
        elseif ( boxInGamma(k) == -1 )
            inGamma_cheb(:,:,j) = false; inGamma_lin(:,:,j) = false;
        end
        % Gamma' locations:
        if ( boxInGamma1(k) == 1 )
            inGamma1_cheb(:,:,j) = true;  inGamma1_lin(:,:,j) = true;
        elseif ( boxInGamma1(k) == -1 )
            inGamma1_cheb(:,:,j) = false; inGamma1_lin(:,:,j) = false;
        end
        % Assemble the point list for the boxes we still don't know.
        % This is so that we may call isinterior() once.
        domk = domains(:,k);
        sclx = diff(domk(1:2));
        scly = diff(domk(3:4));
        xx(:,:,j) = sclx*xx0 + domk(1); xxx(:,:,j) = sclx*xxx0 + domk(1);
        yy(:,:,j) = scly*yy0 + domk(3); yyy(:,:,j) = scly*yyy0 + domk(3);
        j = j + 1;
    end

    inGamma_cheb(:,:,boxesToCheck)   = isinterior(Gamma,   xx(:,:,boxesToCheck),   yy(:,:,boxesToCheck));
    inGamma1_cheb(:,:,boxesToCheck1) = isinterior(rGamma1, xx(:,:,boxesToCheck1),  yy(:,:,boxesToCheck1));
    inStrip_cheb = inGamma_cheb & ~inGamma1_cheb;
    vals = zeros(n_box, n_box, ncheck);
    vals(inGamma_cheb) = f(xx(inGamma_cheb), yy(inGamma_cheb));
    
    if ( strcmpi(errorType, 'vals') )
        inGamma_lin(:,:,boxesToCheck)   = isinterior(Gamma,   xxx(:,:,boxesToCheck),  yyy(:,:,boxesToCheck));
        inGamma1_lin(:,:,boxesToCheck1) = isinterior(rGamma1, xxx(:,:,boxesToCheck1), yyy(:,:,boxesToCheck1));
        inStrip_lin  = inGamma_lin  & ~inGamma1_lin;
        F = zeros(m_box, m_box, ncheck);
        F(inGamma_lin) = f(xxx(inGamma_lin), yyy(inGamma_lin));
    end

    % Now evaluate at the strip points
    j = 1;
    nboxes_old = nboxes;
    for k = id:nboxes_old

        inStrip_cheb_j = inStrip_cheb(:,:,j);
        xx_j = xx(:,:,j);     xx_sj = xx_j(inStrip_cheb_j);
        yy_j = yy(:,:,j);     yy_sj = yy_j(inStrip_cheb_j);
        vals_j = vals(:,:,j); vals_sj = vals_j(inStrip_cheb_j);
        rr_sj  = zeros(size(xx_sj));
        found_cheb = false(size(xx_sj));

        if ( strcmpi(errorType, 'vals') )
            inStrip_lin_j  = inStrip_lin(:,:,j);
            xxx_j = xxx(:,:,j);   xxx_sj = xxx_j(inStrip_lin_j);
            yyy_j = yyy(:,:,j);   yyy_sj = yyy_j(inStrip_lin_j);
            F_j = F(:,:,j);       F_sj = F_j(inStrip_lin_j);
            rrr_sj = zeros(size(xxx_sj));
            found_lin  = false(size(xxx_sj));
        end

        for e = elems{k}
            % Assume the element list is already filtered
            % Filter out the points that are not in the bounding box of this element
            box = [min(dom(e).x(:))-1e-12 max(dom(e).x(:))+1e-12;
                   min(dom(e).y(:))-1e-12 max(dom(e).y(:))+1e-12];         
            inBoundingBox = xx_sj > box(1,1) & xx_sj < box(1,2) & ...
                            yy_sj > box(2,1) & yy_sj < box(2,2);
            inElem = inBoundingBox & ~found_cheb;
            if ( any(inElem) )
                nx = gam_nx_cheb{e}(:,1);
                ny = gam_nx_cheb{e}(:,2);
                [~, rr, found] = getElemCoordinates(xx_sj(inElem), ...
                    yy_sj(inElem), gamx_cheb{e}, gamy_cheb{e}, ...
                    gam1x_cheb{e}, gam1y_cheb{e}, nx, ny, xcheb, vcheb, nthreads);
                found(abs(rr) > 1) = false;
                inElem(inElem) = found;
                rr_sj(inElem) = rr(found);
                found_cheb(inElem) = true;
            end

            if ( strcmpi(errorType, 'vals') )
                % Filter out the points that are not in the bounding box of this element
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
        end

        vals_sj = step(rr_sj) .* vals_sj;
        vals_sj(~found_cheb) = 0;
        vals_j(inStrip_cheb_j) = vals_sj; % This is like: vals = f(xx,yy);
        vals(:,:,j) = vals_j;

        switch lower(errorType)
            case 'vals'
                F_sj = step(rrr_sj) .* F_sj;
                F_sj(~found_lin) = 0;
                F_j(inStrip_lin_j) = F_sj; % This is like: F = f(xxx,yyy)
                F(:,:,j) = F_j;
                G = treefun2.chebvals2refvals(vals_j);
                err = norm(F_j(:) - G(:), inf);
            case 'coeffs'
                coeffs_j = treefun2.vals2coeffs(vals_j);
                Ex = sum(abs(coeffs_j(end-1:end,:)), 'all') / (2*n_box);
                Ey = sum(abs(coeffs_j(:,end-1:end)), 'all') / (2*n_box);
                err = (Ex + Ey) / 2;
            otherwise
                error('Unknown error type.');
        end
        %vscale = max( abs( vals_j(:) ) );
        resolved = ( err < tol );

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

            children(:,k) = [cid1 cid2 cid3 cid4];
            heights(k) = 1;
            nboxes = nboxes + 4;
        end

        j = j + 1;
    end
    
    id = id + ncheck;

    if ( nboxes > maxBoxes )
        warning('Number of boxes has exceeded %d.', maxBoxes);
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
tf = treefun2(domains, levels, heights, ids, parents, children, coeffs, ...
    cols, rows);

% boxes = struct();
% for k = 1:nboxes
%     boxes(k).domain = domains(:,k).';
%     boxes(k).level = levels(k);
%     boxes(k).height = heights(k);
%     boxes(k).id = ids(k);
%     boxes(k).parent = parents(k);
%     boxes(k).children = children(:,k).';
%     boxes(k).coeffs = coeffs{k};
%     boxes(k).col = cols(k);
%     boxes(k).row = rows(k);
% end
% tf = treefun2_old(boxes);

fprintf('   Constructing quadtree .............. %.6fs\n', toc);

end