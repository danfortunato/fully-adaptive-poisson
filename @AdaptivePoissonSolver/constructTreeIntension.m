function [tf, isource] = constructTreeIntension(S, f, opts)

defaults = [];
defaults.debug = true;
defaults.tol = 1e-11;
defaults.boxToStripRatio = 0.5;
defaults.maxBoxes = 300000;

if ( nargin < 3 )
    opts = defaults;
else
    opts = setDefaults(opts, defaults);
end

Gamma  = S.Gamma_re;
Gamma1 = S.Gamma1;
n_box  = S.n_box;
%n_box_alias = 2*n_box; % Need to account for aliasing in the coefficients?
n_box_alias = n_box;
[xx0, yy0] = chebpts2(n_box_alias, n_box_alias, [0 1 0 1]);

% Allocate space for many boxes
maxBoxes = opts.maxBoxes;
domains  = zeros(4, maxBoxes);
levels   = zeros(maxBoxes, 1);
heights  = zeros(maxBoxes, 1);
ids      = zeros(maxBoxes, 1);
parents  = zeros(maxBoxes, 1);
children = zeros(4, maxBoxes);
coeffs   = cell(maxBoxes, 1);
cols     = zeros(maxBoxes, 1, 'uint64');
rows     = zeros(maxBoxes, 1, 'uint64');
needsEval = [];

% Flags to tell where the boxes lie in relation to each curve:
%    1 = fully inside
%   -1 = fully outside
%    0 = neither
boxInGamma  = zeros(maxBoxes, 1);
boxInGamma1 = zeros(maxBoxes, 1);

% Case 1: Fully outside Gamma and Gamma'
% Case 2: Fully inside Gamma and Gamma'
% Case 3: Fully inside Gamma, fully outside Gamma'
% Case 4: Overlapping Gamma, fully outside Gamma'
% Case 5: Overlapping Gamma', fully inside Gamma
% Case 6: Overlapping both Gamma and Gamma'
%
% There are also some impossible cases. Not sure we need to check for
% these:
% - Fully outside Gamma, fully inside Gamma'
% - Overlapping Gamma, fully inside Gamma'
% - Overlapping Gamma', fully outside Gamma

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
cols(1)      = 0;
rows(1)      = 0;

t_interior = 0;
t_near = 0;
t_stripIsResolved = 0;
t_feval = 0;
t_functionIsResolved = 0;
t_split = 0;

% Note: the length changes at each iteration here!
% We try to avoid calling isinterior() on the dense set of leaf points.
% We will proceed in three phases.
id = 1;
while ( id <= nboxes )

    % We will check all the leaves on this level for refinement.
    ilev = id:nboxes;
    nblev = length(ilev);
    dbprint(opts.debug, 'nboxes = %d\n', nboxes);

    %% Phase 1
    %  If our parent is fully inside or outside of both curves then we know
    %  where we are.

    if ( id > 1 )
        boxInGamma(ilev)  = boxInGamma(parents(ilev));
        boxInGamma1(ilev) = boxInGamma1(parents(ilev));
    end

    %% Phase 2
    %  If we weren't able to determine where we are from Phase 1, then we
    %  check where the four corners of our box lie. If they are all on the
    %  same side of the curve and no boundaries lie near the box, then we
    %  know where we are. We use a KD-tree to perform this range query.

    % Collect all the corners of the questionable boxes into an array so we
    % can call isinterior() once, on the questionable ones.
    t = tic;
    ques  = ( boxInGamma(ilev)  == 0 );
    ques1 = ( boxInGamma1(ilev) == 0 );
    corners_x = domains([1 2 2 1], ilev);
    corners_y = domains([3 3 4 4], ilev);
    cornersInGamma  = isinterior(Gamma,  corners_x(:,ques),  corners_y(:,ques));
    cornersInGamma1 = isinterior(Gamma1, corners_x(:,ques1), corners_y(:,ques1));
    t_interior = t_interior + toc(t);

    t = tic;
    j = 1; % Too many indices...
    l = 1;
    l1 = 1;
    for k = id:nboxes
        % We will look for any points in Gamma or Gamma' that lie "near"
        % this box, where "near" means inside of a slightly larger box
        % enclosing this one.
        if ( boxInGamma(k) == 0 || boxInGamma1(k) == 0 )
            % Pad by a multiple of the diagonal length of the box
            pad = 0.1*sqrt(2)*diff(domains(1:2, k));
            box = [domains(1,k)-pad domains(2,k)+pad;
                   domains(3,k)-pad domains(4,k)+pad];
        end
        if ( boxInGamma(k) == 0 )
            % This box's relation to Gamma is still unknown. Check if all
            % the corners lie on the same side and use the KD tree to test
            % for nearness.
            nearGamma = ~isempty(kd.range(box));
            if ( all(cornersInGamma(:,l)) && ~nearGamma )
                boxInGamma(k) = 1;
            elseif ( all(~cornersInGamma(:,l)) && ~nearGamma )
                boxInGamma(k) = -1;
            end
            l = l+1;
        end
        if ( boxInGamma1(k) == 0 )
            % This box's relation to Gamma is still unknown. Check if all
            % the corners lie on the same side and use the KD tree to test
            % for nearness.
            nearGamma1 = ~isempty(kd1.range(box));
            if ( all(cornersInGamma1(:,l1)) && ~nearGamma1 )
                boxInGamma1(k) = 1;
            elseif ( all(~cornersInGamma1(:,l1)) && ~nearGamma1 )
                boxInGamma1(k) = -1;
            end
            l1 = l1+1;
        end
        j = j+1;
    end
    t_near = t_near + toc(t);

    %% We now check to see if the refinement criterion is satisfied in the
    %  strip. Roughly, this criterion says that a box cut by Gamma cannot
    %  be neighbors with a box cut by Gamma'.

    t = tic;
    stripIsResolved = true(nblev, 1);
    big  = boxInGamma(ilev);
    big1 = boxInGamma1(ilev);

    % (1) If the box overlaps the entire strip, then it needs to be further
    %     refined and the strip is not resolved.
    filt = ( big == 0 & big1 == 0 );
    if ( any(filt) )
        stripIsResolved(filt) = false;
    end

    % (2) If the box is fully inside or half overlapping the strip, then we
    %     check if the box size is a fraction of the strip width. If it is,
    %     then the strip is resolved.
    filt = (big == 1 & big1 == -1) | (big == 0 & big1 == -1) | (big == 1 & big1 == 0);
    if ( any(filt) )
        doms = domains(:,ilev(filt));
        boxSize = sqrt(2)*diff(doms(1:2,:)).';
        boxCenter = [0.5*(doms(1,:)+doms(2,:)) ; 0.5*(doms(3,:)+doms(4,:))].';
        idx  = kd.nn(boxCenter);
        idx1 = kd1.nn(boxCenter);
        stripWidth = vecnorm(S.gamxy_re(idx,:) - S.gam1xy(idx1,:), 2, 2);
        fac = zeros(size(big(filt)));
        fac(big(filt) == 0) =   opts.boxToStripRatio;
        fac(big(filt) ~= 0) = 2*opts.boxToStripRatio;
        stripIsResolved(filt) = ( boxSize < fac.*stripWidth );
    end

    t_stripIsResolved = t_stripIsResolved + toc(t);

    %% Now we fill in the function values on the uncut boxes.

    % Evaluate the given function on the uncut interior boxes for which the
    % refinement criterion is satisfied.
    t = tic;
    filt = ( boxInGamma(ilev) == 1 & stripIsResolved );
    nfilt = length(find(filt));
    doms = domains(:,ilev(filt));
    sclx = diff(doms(1:2,:));
    scly = diff(doms(3:4,:));
    xx = sclx.*xx0(:) + doms(1,:); xx = reshape(xx, [n_box_alias n_box_alias nfilt]);
    yy = scly.*yy0(:) + doms(3,:); yy = reshape(yy, [n_box_alias n_box_alias nfilt]);
    if ( isempty(xx) )
        vals = zeros(size(xx));
    else
        vals = f(xx,yy);
    end
    t_feval = t_feval + toc(t);

    j = 1;
    nboxes_old = nboxes;
    for k = id:nboxes_old

        t = tic;
        % Now check if the function is resolved on each box. We can save
        % some work by skipping this step when the strip is not resolved.
        if ( stripIsResolved(k-id+1) )
            if ( boxInGamma(k) == -1 )
                % If the box is fully outside of the domain, it is zero and
                % thus resolved.
                functionIsResolved = true;
                coeffs_j = zeros(n_box);
            elseif ( boxInGamma(k) == 1 )
                % The box does not overlap Gamma, so its function values
                % are all given by the righthand side f.
                vals_j = vals(:,:,j);
                coeffs_j = treefun2.vals2coeffs(vals_j);
                coeffs_j = coeffs_j(1:n_box,1:n_box);
                Ex = sum(abs(coeffs_j(end-1:end,:)), 'all') / (2*n_box);
                Ey = sum(abs(coeffs_j(:,end-1:end)), 'all') / (2*n_box);
                err = (Ex + Ey) / 2;
                vscale = max(abs(vals_j(:)));
                functionIsResolved = ( err < opts.tol * max(vscale, 1) );
                j = j+1;
            else
                % Some of the box's nodes lie outside of the domain.
                % The strip is resolved, so we just need to construct the
                % (discontinuous) function by evaluating the righthand side
                % f at the interior nodes. We do not care if the function
                % is resolved on these boxes.
                functionIsResolved = true;
                needsEval(end+1) = k; %#ok<AGROW>
                coeffs_j = 0;
                % We will evaluate the box coefficients after the loop is finished.
            end
        end
        t_functionIsResolved = t_functionIsResolved + toc(t);

        t = tic;
        if ( stripIsResolved(k-id+1) && functionIsResolved )
            coeffs{k} = coeffs_j;
            heights(k) = 0;
        else
            % Split into four child boxes
            domp = domains(:,k);
            xmid = 0.5*(domp(1)+domp(2));
            ymid = 0.5*(domp(3)+domp(4));

            cid1 = nboxes+1;
            domains(:,cid1) = [domp(1) xmid domp(3) ymid];
            ids(cid1)       = cid1;
            parents(cid1)   = k;
            levels(cid1)    = levels(k)+1;
            heights(cid1)   = 0;
            cols(cid1)      = 2*cols(k);
            rows(cid1)      = 2*rows(k);

            cid2 = nboxes+2;
            domains(:,cid2) = [xmid domp(2) domp(3) ymid];
            ids(cid2)       = cid2;
            parents(cid2)   = k;
            levels(cid2)    = levels(k)+1;
            heights(cid2)   = 0;
            cols(cid2)      = 2*cols(k) + 1;
            rows(cid2)      = 2*rows(k);

            cid3 = nboxes+3;
            domains(:,cid3) = [domp(1) xmid ymid domp(4)];
            ids(cid3)       = cid3;
            parents(cid3)   = k;
            levels(cid3)    = levels(k)+1;
            heights(cid3)   = 0;
            cols(cid3)      = 2*cols(k);
            rows(cid3)      = 2*rows(k) + 1;

            cid4 = nboxes+4;
            domains(:,cid4) = [xmid domp(2) ymid domp(4)];
            ids(cid4)       = cid4;
            parents(cid4)   = k;
            levels(cid4)    = levels(k)+1;
            heights(cid4)   = 0;
            cols(cid4)      = 2*cols(k) + 1;
            rows(cid4)      = 2*rows(k) + 1;

            children(:,k) = [cid1 cid2 cid3 cid4];
            heights(k) = 1;
            nboxes = nboxes + 4;
        end
        t_split = t_split + toc(t);

    end

    id = id + nblev;

    if ( nboxes > 100000 )
        warning('Number of boxes has exceeded %d.', maxBoxes);
        keyboard
    end

end

Timer.log(sprintf('Box interior queries: %gs', t_interior));
Timer.log(sprintf('Box nearness queries: %gs', t_near));
Timer.log(sprintf('Checking strip criteria: %gs', t_stripIsResolved));
Timer.log(sprintf('Checking function resolved: %gs', t_functionIsResolved));
Timer.log(sprintf('Splitting boxes: %gs', t_split));
Timer.log(sprintf('Evaluating uncut boxes: %gs', t_feval));

%% Do a cumulative sum in reverse to correct the heights
for k = nboxes:-1:1
    if ( heights(k) ~= 0 )
        heights(k) = 1 + max(heights(children(:,k)));
    end
end

%% Phase 3
%  Evaluate the function on the cut boxes.
Timer.tic();

% Collect box nodes from the boxes that need evaluation:
xx = zeros(n_box_alias, n_box_alias, length(needsEval));
yy = zeros(n_box_alias, n_box_alias, length(needsEval));
for k = 1:length(needsEval)
    id = needsEval(k);
    xx(:,:,k) = diff(domains(1:2,id))*xx0 + domains(1,id);
    yy(:,:,k) = diff(domains(3:4,id))*yy0 + domains(3,id);
end

% Call isinterior() once:
inGamma = isinterior(Gamma, xx, yy);

% Assign the box coefficients:
vals = zeros(n_box_alias, n_box_alias, length(needsEval));
vals(inGamma) = f(xx(inGamma), yy(inGamma));
for k = 1:length(needsEval)
    id = needsEval(k);
    coeffs{id} = treefun2.vals2coeffs(vals(:,:,k));
    coeffs{id} = coeffs{id}(1:n_box,1:n_box);
end

Timer.toc('Evaluating cut boxes');

%% Convert to treefun2 (this will balance the tree)
Timer.tic();

domains  = domains(:,1:nboxes);
levels   = levels(1:nboxes);
heights  = heights(1:nboxes);
ids      = ids(1:nboxes);
parents  = parents(1:nboxes);
children = children(:,1:nboxes);
coeffs   = coeffs(1:nboxes);
cols     = cols(1:nboxes);
rows     = rows(1:nboxes);

tf = treefun2(domains, levels, heights, ids, parents, children, coeffs, cols, rows);

Timer.toc('Balancing tree');

%% Mark the boxes that are entirely outside of the domain.
%  This is useful when computing the volume potential, as the contribution
%  from these boxes is zero and may be skipped.
Timer.tic();

lid = leaves(tf);
corners_x = tf.domain([1 2 2 1], lid);
corners_y = tf.domain([3 3 4 4], lid);
isource = any(isinterior(Gamma, corners_x, corners_y));

Timer.toc('Marking exterior boxes');

end
