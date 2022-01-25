function z = perturbBreaks(dom, width)
%PERTURBBREAKS   Perturb the breaks in a panelization based on local panel
% size.

if ( nargin < 2 )
    width = 0.3;
end

%nx = normal(chebfun(dom.f, [0 2*pi], 'trig'));
nx = normal(chebfun(dom.f, [0 2*pi]));
nx = nx(:,1)+nx(:,2)*1i;
nx = nx./abs(nx);
normals = nx(dom.breaks(1:end-1)).';

panelsize = sum([dom.w{:}]).';
panelsize = [panelsize(end); panelsize];
%panelsize = (panelsize(1:end-1) + panelsize(2:end)) / 2;
panelsize = sqrt(panelsize(1:end-1).*panelsize(2:end));

z = dom.zbreaks(1:end-1).' - width*panelsize.*normals;

end
