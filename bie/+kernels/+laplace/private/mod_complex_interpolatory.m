function w = mod_complex_interpolatory(tleg, z)

% Helsing & Ojala
N = numel(tleg);
A=ones(N);
for k=2:N
    A(:,k)=tleg.*A(:,k-1);
end
p=complex(zeros(N,1));
p(1)=log((1-z)/(-1-z));
s = -1;
for k=1:N-1
    p(k+1)=z*p(k)+(1-s)/k;
    s = -s; % s=(-1)^k
end
warning('off', 'MATLAB:nearlySingularMatrix')
w = (A.'\p) .* (tleg-z);
warning('on', 'MATLAB:nearlySingularMatrix')

end
