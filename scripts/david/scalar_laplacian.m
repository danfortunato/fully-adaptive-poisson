function luh = scalar_laplacian(CO, AAG, RAG, uh)
    uh_t  = CO.R01 * (uh .* AAG.iks.');
    uh_tt = CO.R12 * (fourier_multiply(uh_t, RAG.inv_psi1) .* AAG.iks.');
    uh_rr = CO.D12 * (fourier_multiply(CO.D01*uh, RAG.psi1));
    luh = fourier_multiply(uh_rr+uh_tt, RAG.inv_psi2);
end

function out = mfft(f)
    [m,n] = size(f);
    n2 = floor(n/2);
    fh = fft(f, [], 2);
    out = zeros(m,n-1);
    out(:,1:n2) = fh(:,1:n2);
    out(:,n2+1:end) = fh(:,n2+2:end);
end

function out = mifft(fh)
    [m,ns] = size(fh);
    n = ns+1;
    n2 = floor(n/2);
    out = zeros(m,n);
    out(:,1:n2) = fh(:,1:n2);
    out(:,n2+2:end) = fh(:,n2+1:end);
    out = ifft(out, [], 2);
end

function out = fourier_multiply(fh, m)
    out = mfft(m.*mifft(fh));
end

