clear all; close all; clc;
addpath('utils')

% Recording analytics
total_tic = tic;

% Parameters
nu = 0.001;
Lx = 20;
Ly = 20;
m = 64;
N = m * m;

% IC-specific parameters
initial_conditions = 4;
spectral_timing_results = zeros(1, initial_conditions);

% Main IC loop
for ic = 1:initial_conditions
    tic;
    
    % Spatial grid setup 
    x2 = linspace(-Lx/2, Lx/2, m+1); x = x2(1:m);
    y2 = linspace(-Ly/2, Ly/2, m+1); y = y2(1:m);
    [X, Y] = meshgrid(x, y);

    % Set initial condition
    switch ic
        case 1 % Single Gaussian
            w = exp(-2*X.^2-Y.^2/20);
        case 2 % Colliding vortices
            w = exp(-2*(X+2).^2-Y.^2/20) + exp(-2*(X-2).^2-Y.^2/20);
        case 3 % Multiple vortices
            w = exp(-2*(X+6).^2-Y.^2/10) + exp(-2*(X-6).^2-Y.^2/10)...
                + exp(-2*X.^2-(Y+6).^2/10) + exp(-2*X.^2-(Y-6).^2/10);
        case 4 % Dipole
            w = exp(-2*(X+2).^2-(Y+6).^2/20) - exp(-2*(X-2).^2-(Y+6).^2/20);
    end

    % Spectral setup
    kx = (2*pi/Lx)*[0:(m/2-1) (-m/2):-1]; kx(1) = -10^(-6);
    ky = (2*pi/Ly)*[0:(m/2-1) (-m/2):-1]; ky(1) = -10^(-6);
    [KX, KY] = meshgrid(kx, ky);
    Kderv = KX.^2 + KY.^2;
    iKX = 1i * KX;
    iKY = 1i * KY;

    % Time integration 
    wt = fft2(w);
    wt2 = reshape(wt, N, 1);

    % Time setup
    tspan = 0:4:32; 
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
    
    % Solve ODEs and store results
    [t, wt2sol] = ode45(@(t, wt2) spc_rhs(t, wt2, Kderv, iKX, iKY, nu, m, m), tspan, wt2, options);

    saved_omega = zeros(m, m, length(tspan));
    figure;
    for j = 1:length(tspan)
        wt2_current = wt2sol(j, :)'; 
        wt = reshape(wt2_current, m, m); 
        saved_omega(:, :, j) = real(ifft2(wt)); 
        subplot(3, 3, j); 
        pcolor(x, y, saved_omega(:, :, j));
        shading interp; colorbar; clim([0 1]);
        title(['t = ', num2str(tspan(j))], 'FontSize', 14);
    end

    % Store timing results
    time_taken = toc;
    spectral_timing_results(ic) = time_taken;

    % Process results for IC1 only
    if ic == 1

        % RMSE calculation
        load('spectral_highres.mat');
        highres_size = size(saved_omega_spectral_highres, 1);
        block_size = highres_size/m;
        
        saved_omega_highres_downsampled = zeros(m, m, length(tspan));
        for k = 1:length(tspan)
            for i = 1:m
                for j = 1:m
                    row_start = (i-1)*block_size + 1;
                    row_end = i*block_size;
                    col_start = (j-1)*block_size + 1;
                    col_end = j*block_size;
                    block = real(saved_omega_spectral_highres(row_start:row_end, col_start:col_end, k));
                    saved_omega_highres_downsampled(i,j,k) = mean(block(:), 'omitnan');
                end
            end
        end

        rmse_spectral = zeros(1, length(tspan));
        for k = 1:length(tspan)
            diff = saved_omega(:, :, k) - saved_omega_highres_downsampled(:, :, k);
            rmse_spectral(k) = sqrt(mean(diff(:).^2));
        end
        
        % Save IC1 results
        save(fullfile('data', sprintf('spectral_ic1_%dx%d.mat', m, m)), 'rmse_spectral');
    end
end

% Save aggregate timing
save(fullfile('data', 'spectral_timing_results.mat'), 'spectral_timing_results');
fprintf('Total simulation time: %.2f seconds\n', toc(total_tic));