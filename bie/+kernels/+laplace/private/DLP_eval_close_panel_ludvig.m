function u = DLP_eval_close_panel_ludvig(source, target, sigma, side)
%DLP_EVAL_CLOSE_PANEL   Evaluate double-layer potential at target
%points using panel close evaluation.

u = FMM_eval(source, target, 'dipstr', sigma);

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
            if ( ~isempty(source.nearby_schwarz) && source.nearby_schwarz(k) )
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
                %udirect = sum( imag(den.*source.cw{k}./(source.z{k}-zt)) ) / (2*pi);
                pa = struct();
                pa.x  = source.z{k};
                pa.xp = source.dz{k};
                pa.nx = source.normal{k}(:,1) + 1i*source.normal{k}(:,2);
                pa.sp = source.speed{k};
                pa.w  = source.w{k};
                pa.wxp = source.cw{k};
                D_direct = LapDLP(struct('x', zt), pa);
                udirect = -D_direct * den;

                %wreal = mod_complex_interpolatory(t01, t) * source.w{k}(1)/w01(1) / source.speed{k}(1);
                %wreal = real_quad_2d(t01, t) .* source.w{k}./w01 ./ source.speed{k};
                wreal = real_quad_2d(t01, t) .* source.w{k}(1)/w01(1) / source.speed{k}(1);
                ureal = imag( sum(den.*source.dz{k}./(source.z{k}-zt) .* wreal) ) / (2*pi);

                % Apply
                u(i) = u(i) + udirect - ureal;
            end
        end
    end
end

%u(invalid) = nan;

end
