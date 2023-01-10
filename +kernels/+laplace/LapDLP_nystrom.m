function D = LapDLP_nystrom(source)

n = source.N;
np = source.np;
z = cell2mat(source.z);
t = struct();
t.x = z;

D = zeros(n*np);

for k = 1:np
    % Package source panel into structs for close evaluation routine
    pa = struct();
    pa.x  = source.z{k};
    pa.xp = source.dz{k};
    pa.nx = source.normal{k}(:,1) + 1i*source.normal{k}(:,2);
    pa.sp = source.speed{k};
    pa.w  = source.w{k};
    pa.wxp = source.cw{k};

    selfIdx = (1:n)+(k-1)*n;
    D(:,selfIdx) = D(:,selfIdx) + LapDLP(t, pa);
end

kappa = cell2mat(source.curvature);
w = cell2mat(source.w);
D(1:n*np+1:end) = -kappa/(4*pi).*w;

end
