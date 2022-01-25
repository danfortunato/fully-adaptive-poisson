function [X_sm, Y_sm, h] = smoother2D(XY, n_skel, n_curve)
%SMOOTHER2D   A 2D version of the multiscale surface smoother.

% Controls the base smoothing size
lambda = 0.2;

% Wrap the points to close the curve
XY = [XY; XY(1,:)];

[x0, y0, sigma0, Nx_vert, Ny_vert] = centers_and_values(XY, lambda);

[X_sk, Y_sk, W_sk, Nx_sk, Ny_sk] = setup_integration_nodes(XY, n_skel);

[Xb_sm, Yb_sm, psNx_sm, psNy_sm] = setup_base_points(XY, n_curve, Nx_vert, Ny_vert);

[X_sm, Y_sm, h] = mi_interp_newton(Xb_sm, Yb_sm, psNx_sm, psNy_sm, X_sk, Y_sk, W_sk, Nx_sk, Ny_sk, x0, y0, sigma0);

end

function [x0, y0, sigma, Nx_vert, Ny_vert] = centers_and_values(XY, lambda)

    m = size(XY,1) - 1;
    PP = (XY(1:m,:) + XY(2:m+1,:)) / 2;
    x0 = PP(:,1);
    y0 = PP(:,2);
    d = sqrt((XY(1:m,1)-XY(2:m+1,1)).^2 + (XY(1:m,2)-XY(2:m+1,2)).^2);
    %fprintf('mid(d) = %g\n', min(d))
    sigma = d*lambda;

    P1 = XY(1:m,:);
    P2 = XY(2:m+1,:);
    Pdiff = P2 - P1;
    Norm_x =  Pdiff(:,2);
    Norm_y = -Pdiff(:,1);
    Norm_x = [Norm_x(end); Norm_x];
    Norm_y = [Norm_y(end); Norm_y];
    d_Norm = sqrt(Norm_x.^2 + Norm_y.^2);
    Norm_x = Norm_x ./ d_Norm;
    Norm_y = Norm_y ./ d_Norm;
    Nx_vert = (Norm_x(1:m) + Norm_x(2:m+1)) / 2;
    Ny_vert = (Norm_y(1:m) + Norm_y(2:m+1)) / 2;

end

function [X_sk, Y_sk, W_sk, Nx_sk, Ny_sk] = setup_integration_nodes(XY, n)

    %[t, w] = legpts(n, [0 1]);
    [t, w] = chebpts(n, [0 1]);
    t = t.';

    m = size(XY,1) - 1;
    P1 = XY(1:m,:);
    P2 = XY(2:m+1,:);
    Pdiff = P2 - P1;
    Norm_x =  Pdiff(:,2);
    Norm_y = -Pdiff(:,1);
    dl = sqrt((XY(1:m,1)-XY(2:m+1,1)).^2 + (XY(1:m,2)-XY(2:m+1,2)).^2);
    Nx_sk = (Norm_x ./ dl) * ones(size(t));
    Ny_sk = (Norm_y ./ dl) * ones(size(t));
    X_sk = P1(:,1) + (P2(:,1)-P1(:,1))*t;
    Y_sk = P1(:,2) + (P2(:,2)-P1(:,2))*t;
    W_sk = dl*w;

    X_sk  = reshape(X_sk',  [], 1);
    Y_sk  = reshape(Y_sk',  [], 1);
    W_sk  = reshape(W_sk',  [], 1);
    Nx_sk = reshape(Nx_sk', [], 1);
    Ny_sk = reshape(Ny_sk', [], 1);

end

function [Xb_sm, Yb_sm, psNx_sm, psNy_sm] = setup_base_points(XY, n, Nx_sk, Ny_sk)

    %t = legpts(n, [0 1]);
    t = chebpts(n, [0 1]);
    t = t.';

    m = size(XY,1) - 1;
    P1 = XY(1:m,:);
    P2 = XY(2:m+1,:);
    Xb_sm = P1(:,1) + (P2(:,1)-P1(:,1))*t;
    Yb_sm = P1(:,2) + (P2(:,2)-P1(:,2))*t;
    Xb_sm = reshape(Xb_sm', [], 1);
    Yb_sm = reshape(Yb_sm', [], 1);

    Nx_sk = [Nx_sk; Nx_sk(1)];
    Ny_sk = [Ny_sk; Ny_sk(1)];

    psNx_sm = Nx_sk(1:m)*(1-t) + Nx_sk(2:m+1)*t;
    psNy_sm = Ny_sk(1:m)*(1-t) + Ny_sk(2:m+1)*t;

    psNx_sm = reshape(psNx_sm', [], 1);
    psNy_sm = reshape(psNy_sm', [], 1);

end

function [sigma, dx_sigma, dy_sigma] = sigma_eval(X, Y, x0, y0, sigma0)

    di = (x0-X').^2 + (y0-Y').^2;
    alpha = 1./(5*max(sigma0).^2)/2;
    %fprintf('alpha = %g\n', alpha)
    DDExp = exp(-alpha.*di);
    F = sigma0.' * DDExp;
    D = sum(DDExp);
    %fprintf('min(F) = %g, max(D) = %g\n', min(F), max(D))
    sigma = (F ./ D)';
    dxDDExp = -alpha*2*(x0-X') .* DDExp;
    dyDDExp = -alpha*2*(y0-Y') .* DDExp;
    dFdx = sigma0.' * dxDDExp;
    dDdx = sum(dxDDExp);
    dFdy = sigma0.' * dyDDExp;
    dDdy = sum(dyDDExp);
    dx_sigma = (-(dFdx.*D - F.*dDdx) ./ D.^2)';
    dy_sigma = (-(dFdy.*D - F.*dDdy) ./ D.^2)';

end

function [F, gradF_x, gradF_y] = eval_F(X, Y, X_sk, Y_sk, W_sk, Nx_sk, Ny_sk, x0, y0, sigma0)

    [sigma, dx_sigma, dy_sigma] = sigma_eval(X, Y, x0, y0, sigma0);
    %fprintf('min(sigma0) = %g, min(sigma) = %g\n', min(sigma0), min(sigma))

    d2 = (X_sk-X').^2 + (Y_sk-Y').^2;
    [DDD, DDDp, DDDsigma] = mi_phi_der(d2, sigma);

    F = (Nx_sk.*W_sk)'*(DDD.*(X_sk-X')) + (Ny_sk.*W_sk)'*(DDD.*(Y_sk-Y'));
    F = F';
    a = length(dx_sigma);
    TTT_dx_sigma = sparse(1:a, 1:a, dx_sigma);
    TTT_dy_sigma = sparse(1:a, 1:a, dy_sigma);

    gradF_x_aux = -DDDp.*(X_sk-X') + DDDsigma*TTT_dx_sigma;
    gradF_x_1 = (Nx_sk.*W_sk)'*(gradF_x_aux.*(X_sk-X')) + (Ny_sk.*W_sk)'*(gradF_x_aux.*(Y_sk-Y'));
    gradF_x_2 = -(Nx_sk.*W_sk)'*DDD;
    gradF_x = (gradF_x_1 + gradF_x_2)';

    gradF_y_aux = -DDDp.*(Y_sk-Y') + DDDsigma*TTT_dy_sigma;
    gradF_y_1 = (Nx_sk.*W_sk)'*(gradF_y_aux.*(X_sk-X')) + (Ny_sk.*W_sk)'*(gradF_y_aux.*(Y_sk-Y'));
    gradF_y_2 = -(Ny_sk.*W_sk)'*DDD;
    gradF_y = (gradF_y_1 + gradF_y_2)';

    F = F - 0.5;

end

function [F, Fp, Fsigma] = mi_phi_der(r2, sigma)
    %This is (1/r)*(d Phi/dr)
    a = length(sigma);
    TTT = sparse(1:a, 1:a, 1./(2*sigma.^2));
    TTT3=sparse(1:a, 1:a, 1./(sigma.^3));
    Mexp = exp(-r2*TTT);
    x = r2*TTT;
    y = expm1ox(x);
    F = -y*TTT/(2*pi);

    y2 = expm1pxox2(x);
    Fp = y2*(TTT*TTT)/pi;

    Fsigma = Mexp*TTT3/(2*pi);
    Fsigma = -Fsigma;

end

function y = expm1ox(x)
    y = zeros(size(x));
    idx1 = find(abs(x) < 0.05);
    idx2 = find(abs(x) >= 0.05);
    y(idx2) = (exp(-x(idx2)) - 1) ./ x(idx2);
    y(idx1) = -x(idx1).^10/39916800 + x(idx1).^9/3628800 - x(idx1).^8/362880 + x(idx1).^7/40320 - x(idx1).^6/5040 + x(idx1).^5/720 - x(idx1).^4/120 + x(idx1).^3/24 - x(idx1).^2/6 + x(idx1)/2 - 1;
end

function y = expm1pxox2(x)
    y = zeros(size(x));
    idx1 = find(abs(x) < 0.1);
    idx2 = find(abs(x) >= 0.1);
    y(idx2) = (exp(-x(idx2)) - 1 + x(idx2).*exp(-x(idx2))) ./ x(idx2).^2;
    y(idx1) = x(idx1).^9/3991680 - x(idx1).^8/403200 + x(idx1).^7/45360 - x(idx1).^6/5760 + x(idx1).^5/840 - x(idx1).^4/144 + x(idx1).^3/30 - x(idx1).^2/8 + x(idx1)/3 - .5;
end

function [F, Fp] = mifyfp(h, Xb_sm, Yb_sm, psNx_sm, psNy_sm, X_sk, Y_sk, W_sk, Nx_sk, Ny_sk, x0, y0, sigma0, idx)
    X = Xb_sm(idx) + psNx_sm(idx).*h(idx);
    Y = Yb_sm(idx) + psNy_sm(idx).*h(idx);
    [F, gradF_x, gradF_y] = eval_F(X, Y, X_sk, Y_sk, W_sk, Nx_sk, Ny_sk, x0, y0, sigma0);
    Fp = gradF_x.*psNx_sm(idx) + gradF_y.*psNy_sm(idx);
end

function [X_sm, Y_sm, h] = mi_interp_newton(Xb_sm, Yb_sm, psNx_sm, psNy_sm, X_sk, Y_sk, W_sk, Nx_sk, Ny_sk, x0, y0, sigma0)

    tol = 1e-10;
    maxiter = 30;

    h  = zeros(size(Xb_sm));
    f  = zeros(size(Xb_sm));
    fp = zeros(size(Xb_sm));

    err = tol+1;
    niter = 0;
    idx = (1:length(h))';

    while ( niter <= maxiter && err > tol )
        [f(idx), fp(idx)] = mifyfp(h, Xb_sm, Yb_sm, psNx_sm, psNy_sm, X_sk, Y_sk, W_sk, Nx_sk, Ny_sk, x0, y0, sigma0, idx);
        delt = f ./ fp;
        h = h - delt;
        errvec = abs(delt);
        idx = find(errvec > tol);
        err = max(max(errvec));
        niter = niter+1;
    end

    X_sm = Xb_sm + psNx_sm.*h;
    Y_sm = Yb_sm + psNy_sm.*h;

end
