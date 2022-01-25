Gamma = Boundary.wavepacket(n, 'quadrature', 'panel');
z = cell(Gamma.np, 1);
width = 0.01;

for k = 1:Gamma.np
    nx = Gamma.normal{k};
    nx = nx(:,1)+nx(:,2)*1i;
    nx = nx./abs(nx);
    z{k} = Gamma.z{k} - width*nx;
end

Gamma1 = Boundary(z);

plot(Gamma), hold on
plot(Gamma.zbreaks, 'k.', 'MarkerSize', 25)
plot(Gamma1, 'r'), hold on
plot(Gamma1.zbreaks, 'r.', 'MarkerSize', 25)
axis equal tight off
shg

%%
plot(Gamma), hold on
plot(Gamma.zbreaks, 'k.', 'MarkerSize', 25)
plot(Gamma1, 'r'), hold on
plot(Gamma1.zbreaks, 'r.', 'MarkerSize', 25)
axis equal tight off
shg