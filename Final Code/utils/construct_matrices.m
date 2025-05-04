function [A, B, C] = construct_matrices(m, delta)
    % Construct FDM matrices
    N = m^2; % Total number of grid points
    Nx = m;  % Number of grid points in x-direction
    
    % Create vectors for matrix construction
    e0 = zeros(N, 1); % Vector of zeros
    e1 = ones(N, 1);  % Vector of ones
    e2 = e1;          % Copy the ones vector
    e4 = e0;          % Copy the zeros vector
    
    for j = 1:Nx
        e2(Nx*j) = 0; % Overwrite every Nx-th value with zero
        e4(Nx*j) = 1; % Overwrite every Nx-th value with one
    end
    
    e3 = circshift(e2, 1); % Shift e2 to the right
    e5 = circshift(e4, 1); % Shift e4 to the right
    
    % Create matrix A (Laplacian)
    A = spdiags([e1 e1 e5 e2 -4*e1 e3 e4 e1 e1], ...
        [-(N-Nx) -Nx -Nx+1 -1 0 1 Nx-1 Nx (N-Nx)], N, N);
    A(1, 1) = 2; % Adjust boundary condition
    A = A / (delta^2); % Scale by grid spacing
    
    % Create matrix B (first-order derivative in x)
    adds = ((1/2)/delta) * ones(N, 1);
    subs = -adds;
    B = spdiags([adds adds subs subs], [Nx -(N-Nx) (N-Nx) -Nx], N, N);
    
    % Create matrix C (first-order derivative in y)
    C = spdiags([e5 -e2 e3 -e4], [-(Nx-1) -1 1 (Nx-1)], N, N);
    C = ((1/2)/delta) .* C;
end