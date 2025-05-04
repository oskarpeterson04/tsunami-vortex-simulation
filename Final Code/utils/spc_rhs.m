function rhs = spc_rhs(~, wt2, Kderv, iKX, iKY, nu, nx, ny)
    wt = reshape(wt2, nx, ny); % Reshape to 2D
    w = real(ifft2(wt)); % Physical space vorticity
    psi = real(ifft2(-wt ./ Kderv)); % Stream function in physical space

    % Compute derivatives in Fourier space
    wy = real(ifft2(iKY .* wt));
    wx = real(ifft2(iKX .* wt));
    psiy = real(ifft2(iKY .* fft2(psi)));
    psix = real(ifft2(iKX .* fft2(psi)));

    % Compute RHS
    rhs = -nu * Kderv .* wt - fft2(psix .* wy) + fft2(psiy .* wx);
    rhs = reshape(rhs, nx * ny, 1); % Reshape to column vector
end
