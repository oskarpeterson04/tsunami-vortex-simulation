clear all; close all; clc;
addpath('utils')

% Recording analytics
total_tic = tic;

% Parameters
nu = 0.001; 
m = 64; 
Lx = 20;
Ly = 20;
delta = 20/m;
N = m*m;

% IC-specific parameters
dt_values = [0.05, 0.05, 0.05, 0.05];
initial_conditions = 4;
fdm_timing_results = zeros(1, initial_conditions);

% Main IC loop
for ic = 1:initial_conditions
    tic; % Start timer for this IC
    dt = dt_values(ic);
    tspan = 0:dt:32;
    save_times = 0:4:32;
    save_indices = round(save_times / dt) + 1;

    saved_omega = zeros(m, m, length(save_indices));
    save_count = 1; 

    % Spatial grid
    x2 = linspace(-Lx/2, Lx/2, m + 1); 
    x = x2(1:m);
    y2 = linspace(-Ly/2, Ly/2, m + 1); 
    y = y2(1:m);
    [X,Y] = meshgrid(x,y);

    % Initial condition
    switch ic
        case 1 % Single Gaussian
            w0 = exp(-2*X.^2-Y.^2/20);
        case 2 % Colliding vortices
            w0 = exp(-2*(X+2).^2-Y.^2/20) + exp(-2*(X-2).^2-Y.^2/20);
        case 3 % Multiple vortices
            w0 = exp(-2*(X+6).^2-Y.^2/10) + exp(-2*(X-6).^2-Y.^2/10)...
                + exp(-2*X.^2-(Y+6).^2/10) + exp(-2*X.^2-(Y-6).^2/10);
        case 4 % Dipole
            w0 = exp(-2*(X+2).^2-(Y+6).^2/20) - exp(-2*(X-2).^2-(Y+6).^2/20);
    end

    % Construct matrices 
    [A, B, C] = construct_matrices(m, delta);

    % Time integration
    w_vec = reshape(w0, [N 1]);
    [Time,Omega] = ode45('fdm_rhs',tspan,w_vec,[],A,B,C);

    % Plotting
    figure;
    plot_index = 1;
    for j = 1:length(tspan)
        if ismember(j, save_indices)
            subplot(3, 3, plot_index); % 3x3 grid
            pcolor(x, y, abs(reshape(Omega(j, :), m, m))); % Reshape and plot
            shading interp;
            colorbar;
            clim([0 1]);
            title(['t = ', num2str((j - 1) * dt)], 'FontSize', 14);
            plot_index = plot_index + 1; 
            save_count = save_count + 1;
        end
    end
    
    % Store timing results
    time_taken = toc;
    fdm_timing_results(ic) = time_taken;

    % Process results for IC1 only
    if ic == 1
        % Extract saved results
        for k = 1:length(save_indices)
            saved_omega(:,:,k) = reshape(Omega(save_indices(k), :), m, m);
        end

        % Calculate RMSE
        load('spectral_highres.mat');
        highres_size = size(saved_omega_spectral_highres, 1);
        block_size = highres_size/m;
        
        saved_omega_highres_downsampled = zeros(m, m, length(save_indices));
        for k = 1:length(save_indices)
            for i = 1:m
                for j = 1:m
                    % Calculate block indices
                    row_start = (i-1)*block_size + 1;
                    row_end = i*block_size;
                    col_start = (j-1)*block_size + 1;
                    col_end = j*block_size;
        
                    % Extract and average block
                    block = saved_omega_spectral_highres(row_start:row_end, col_start:col_end, k);
                    saved_omega_highres_downsampled(i,j,k) = mean(block(:), 'omitnan');
                end
            end
        end

        rmse_fdm = zeros(1, length(save_indices));
        for k = 1:length(save_indices)
            diff = saved_omega(:,:,k) - saved_omega_highres_downsampled(:,:,k);
            rmse_fdm(k) = sqrt(mean(diff(:).^2));
        end
        
        % Save IC1 results
        save(fullfile('data', sprintf('fdm_ic1_%dx%d.mat', m, m)), 'rmse_fdm');
    end
end

% Save aggregate timing
save(fullfile('data', 'fdm_timing_results.mat'), 'fdm_timing_results');
fprintf('Total simulation time: %.2f seconds\n', toc(total_tic));
