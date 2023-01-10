function D_corr = DLP_close_panel_mat(source, side)

% Closeness factor
closepan = 1.2;

n = source.N;
np = source.np;
%D_corr = sparse(n*np, n*np);
D_corr = spalloc(n*np, n*np, 3*n*np);

z = cell2mat(source.z);
kd = KDTree([real(z) imag(z)]);

for k = 1:np

    % Package source panel into struct for close evaluation routine
    a = source.zbreaks(k);
    b = source.zbreaks(k+1);
    pa = struct();
    pa.x  = source.z{k};
    pa.xp = source.dz{k};
    pa.nx = source.normal{k}(:,1) + 1i*source.normal{k}(:,2);
    pa.sp = source.speed{k};
    pa.w  = source.w{k};
    pa.wxp = source.cw{k};

    % Find nearby panels using kd-tree
    % nodes = [real(a) imag(a) ; source.x{k} ; real(b) imag(b)];
    % xymin = min(nodes);
    % xymax = max(nodes);
    % box = [xymin(1) xymax(1) xymin(2) xymax(2)];
    % closeNodeIdx = kd.range(box);
    center = (a+b)/2;
    center = [real(center) imag(center)];
    radius = closepan * abs(b-a);
    closeIdx = kd.ball(center, radius);

    % Don't correct self
    selfIdx = (1:n)+(k-1)*n;
    closeIdx = setdiff(closeIdx, selfIdx);

    % For each nearby panel, add the close evaluation matrix and subtract 
    % the direct evaluation matrix. This yields the correction matrix.
    t = struct();
    t.x = z(closeIdx);
    D_close  = LapDLP_closepanel(t, pa, a, b, side);
    D_direct = LapDLP(t, pa);
    D_corr(closeIdx,selfIdx) = D_corr(closeIdx,selfIdx) + D_close - D_direct;

%     closePanelIdx = unique(ceil(closeIdx / n));
%     for j = closePanelIdx(:).'
% 
%         if ( j == k ), continue, end
% 
%         % For each panel j near to panel k, compute the close evaluation
%         % matrix and the direct evaluation matrix
%         t = struct();
%         t.x = source.z{j};
%         D_close  = LapDLP_closepanel(t, pa, a, b, side);
%         D_direct = LapDLP(t, pa);
% 
%         % Add and subtract them to get the correction matrix
%         bj = (1:n)+(j-1)*n;
%         bk = (1:n)+(k-1)*n;
%         D_corr(bj,bk) = D_corr(bj,bk) + D_close - D_direct;
%     end

end

end
