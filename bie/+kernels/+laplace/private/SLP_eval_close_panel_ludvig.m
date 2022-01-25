function u = SLP_eval_close_panel_ludvig(source, target, sigma, side)
%SLP_EVAL_CLOSE_PANEL   Evaluate single-layer potential at target
%points using panel close evaluation.

u = FMM_eval(source, target, 'charge', sigma);

%% Correct for close panels
marker = FMM_eval(source, target, 'dipstr', ones(size(sigma)));
if ( side == 'i' )
    invalid = abs(marker) < 1e-8;
else
    invalid = abs(marker+1) < 1e-8;
end
R = 3;

% Precompute mappings
L = util.legendre_matrix(source.N);
zhat = complex(zeros(source.N, source.np));
for k = 1:source.np
    za = source.zbreaks(k);
    zb = source.zbreaks(k+1);
    zrs = util.rotate_and_scale(za, zb, source.z{k});
    zhat(:,k) = L*zrs;
end
[t01, w01] = util.gauss(source.N);
% Matrix used for computing weights
A = fliplr(vander(t01));

Z = target(:,1) + target(:,2)*1i;
zz = cat(1, source.z{:});
for i = 1:numel(Z)
    if ( invalid(i) ), continue, end
    zt = Z(i);
    % Find closest point
    [~, imin] = min(abs(zt-zz));
    % Translate to panel idx
    nearest_panel_idx = ceil(imin / source.N);
    % Check vs all panels
    for k = 1:source.np

        h = sum( abs(source.cw{k}) ); % Panel length
        d = min(abs(zt-source.z{k}));

        if ( d < h ) % Coarse filter
            % We seem to be close: Find root
            za = source.zbreaks(k);
            zb = source.zbreaks(k+1);
            zr = util.rotate_and_scale(za, zb, zt); % Use pre-rotation
            % The following part (rootfinding) is the
            % bottleneck for this Matlab code
            if ( source.nearby_schwarz(k) )
                % Matrix:
                coeffs = zhat(:,k);
                coeffs(1) = coeffs(1)-zr;
                t_all = legendre.roots(coeffs);
                t = t_all(1);
            else
                % Newton:
                t0 = zr;
                t = util.solve_legendre(zhat(:,k), zr, t0, 20, 1e-13);
            end
            if ( ((side=='i' && imag(t)<0) || (side=='e' && imag(t)>0)) && k == nearest_panel_idx )
                % Deemed outside by nearest panel
                invalid(i) = 1;
                continue
            end
            if ( util.bernstein_rad(t) < R )
                % We're within critical Bernstein
                den = sigma((1:source.N) + source.N*(k-1));

                % Direct interaction for subtraction
                %%udirect = sum( imag(den.*source.cw{k}./(source.z{k}-zt)) ) / (2*pi);
                %udirect = sum( imag(den.*source.cw{k}.*log((source.z{k}-zt))) ) / (2*pi);
                %udirect = sum( (den.*log(abs(zt-source.z{k}))) .* abs(source.dz{k}) .* w01) / (2*pi);
                %udirect = sum( (den.*log(abs(zt-source.z{k}))) .* source.w{k}) / (2*pi);
                pa = struct();
                pa.x  = source.z{k};
                pa.xp = source.dz{k};
                pa.nx = source.normal{k}(:,1) + 1i*source.normal{k}(:,2);
                pa.sp = source.speed{k};
                pa.w  = source.w{k};
                pa.wxp = source.cw{k};
                %S_direct = LapSLP(struct('x', zt), pa);
                %udirect = -S_direct * den;

                %%wreal = mod_complex_interpolatory(t01, t) * source.w{k}(1)/w01(1) / source.speed{k}(1);
                %wreal = real_quad_2d(t01, t) * source.w{k}(1)/w01(1) / source.speed{k}(1);
                %ureal = imag( sum(den.*source.dz{k}./(source.z{k}-zt) .* wreal) ) / (2*pi);

                % Compute complex integrals for t
                [pL, ~] = complex_integrals(t, source.N);
                % Compute weights
                wL = A.' \ real(pL); % Can use LU here since A is constant
                wLcorr = wL - log(abs(t-t01)) .* w01; % Log correction weights
                %wLcorr = wL - log(abs(t-t01)) .* source.w{k} ./ source.speed{k}; % Log correction weights
                %wLcorr = wLcorr .* source.w{k}(1) ./ w01(1) ./ source.speed{k}(1);
                ureal = sum( den .* abs(source.dz{k}) .* wLcorr) / (2*pi);
                %ureal = sum( imag(den.*source.cw{k}.*wLcorr) ) / (2*pi);

                %wreal = wL - log(abs(t-t01)) .* w01 * source.w{k}(1)/w01(1) / source.speed{k}(1);
                %ureal = sum(den.*abs(source.dz{k}).*log(abs(source.z{k}-zt)) .* wreal) / (2*pi);

                %%wreal = wL * source.w{k}(1) / w01(1) / source.speed{k}(1);
                %%wreal = real_quad_2d(t01, t) .* source.w{k}./w01 ./ source.speed{k};
                %wreal = wL .* source.w{k} ./ w01 ./ source.speed{k};
                %ureal = imag( sum(den.*source.dz{k}.*log(abs(zt-source.z{k})) .* wreal) ) / (2*pi);

                % Apply
                u(i) = u(i) - ureal;
            end
        end
    end
end

%u(invalid) = nan;

end
