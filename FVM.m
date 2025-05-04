% Inspired by the work of Karthik Velakur
% Adapted and extended for multiple initial conditions and timing analysis

clear all; close all; clc;
addpath('utils')

% Recording analytics
total_tic = tic;

% Grid setup
m = 64; % Grid resolution
delta = 20/m; % Grid spacing
x = linspace(-10, 10, m); 
y = linspace(-10, 10, m);
[X, Y] = meshgrid(x, y);

% Common parameters
save_times = 0:4:32; 
nu = 0.001;

% IC-specific time steps
dt_values = [0.05, 0.02, 0.01, 0.02]; 

% Initialize storage
initial_conditions = 4;
fvm_timing_results = zeros(1, initial_conditions);

% Main IC loop
for ic = 1:initial_conditions
    tic; % Start timer for this IC
    
    % Set IC-specific time step
    dt = dt_values(ic);
    
    % Recalculate time parameters for this IC
    tspan = 0:dt:32;
    save_indices = round(save_times / dt) + 1;
    
    % Set initial condition (same as before)
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

    % Construct FDM matrices (same for all ICs)
    N = m^2;
    Nx = m;
    [A, B, C] = construct_matrices(m, delta);

    % Initialize variables
    figure;
    plot_index = 1; 
    u_max = 0;
    v_max = 0;
    saved_omega = zeros(m, m, length(save_indices)); 
    save_count = 1; 

    % Time-stepping loop
    for j = 1:length(tspan)
        omega_old = w;

        % Solve Poisson equation for stream function
        psi = reshape(A \ omega_old(:), m, m);
    
        % Compute derivatives of psi
        dpsidx = reshape(B * psi(:), m, m); % x-derivative
        dpsidy = reshape(C * psi(:), m, m); % y-derivative
    
        % Velocity components
        u = -dpsidy; % Cell-centered u
        v = dpsidx;  % Cell-centered v
    
        % Update maximum velocities
        current_u_max = max(abs(u(:)));
        current_v_max = max(abs(v(:)));
        u_max = max(u_max, current_u_max);
        v_max = max(v_max, current_v_max);
    
        % %--------------------------------------------------------------
        % % Central Difference Flux 
        % %--------------------------------------------------------------
        % % Flux in x-direction
        % flux_x = u .* 0.5 .* (omega_old + circshift(omega_old, [0 -1]));
        % 
        % % Flux in y-direction
        % flux_y = v .* 0.5 .* (omega_old + circshift(omega_old, [-1 0]));
        % %--------------------------------------------------------------
    
        %--------------------------------------------------------------
        % QUICK Scheme for Flux Reconstruction with Face Velocities
        %--------------------------------------------------------------
        % Compute face-centered velocities
        u_face_east = 0.5 * (u + circshift(u, [0 -1]));
        v_face_north = 0.5 * (v + circshift(v, [-1 0]));
        
        % QUICK fluxes at east faces
        positive_east = (u_face_east > 0);
        phi_east_pos = omega_old + (1/8)*(-circshift(omega_old, [0 1]) - 2*omega_old + 3*circshift(omega_old, [0 -1]));
        phi_east_neg = circshift(omega_old, [0 -1]) + (1/8)*(3*omega_old - 2*circshift(omega_old, [0 -1]) - circshift(omega_old, [0 -2]));
        flux_east = u_face_east .* (positive_east .* phi_east_pos + ~positive_east .* phi_east_neg);
        
        % QUICK fluxes at north faces
        positive_north = (v_face_north > 0);
        phi_north_pos = omega_old + (1/8)*(-circshift(omega_old, [1 0]) - 2*omega_old + 3*circshift(omega_old, [-1 0]));
        phi_north_neg = circshift(omega_old, [-1 0]) + (1/8)*(3*omega_old - 2*circshift(omega_old, [-1 0]) - circshift(omega_old, [-2 0]));
        flux_north = v_face_north .* (positive_north .* phi_north_pos + ~positive_north .* phi_north_neg);
        
        % Divergence of fluxes
        flux_west = circshift(flux_east, [0 1]);
        flux_south = circshift(flux_north, [1 0]);
        div_x = (flux_east - flux_west) / delta;
        div_y = (flux_north - flux_south) / delta;
        
        % Update vorticity with diffusion
        diffusion_term = nu * reshape(A * omega_old(:), m, m); % Viscous diffusion
        w = omega_old - dt * (div_x + div_y) + dt * diffusion_term;
        
        % Store results at save points
        if ismember(j, save_indices)
            saved_omega(:, :, save_count) = w; % Store current vorticity
            subplot(3, 3, plot_index); 
            pcolor(x, y, abs(w)); 
            shading interp;
            colorbar;
            clim([0 1]);
            title(['t = ', num2str(tspan(j))], 'FontSize', 14);
            plot_index = plot_index + 1; 
            save_count = save_count + 1;
        end
    end

    % Store timing results
    time_taken = toc; % Capture time for this IC
    fvm_timing_results(ic) = time_taken;
    
    % Calculate and save RMSE only for IC1
    if ic == 1
        % Downsample high-res reference data
        load('spectral_highres.mat');
        highres_size = size(saved_omega_spectral_highres, 1);
        block_size = highres_size/m;
        
        saved_omega_highres_downsampled = zeros(m, m, length(save_indices));
        for k = 1:length(save_indices)
            for i = 1:m
                for j = 1:m
                    row_start = (i-1)*block_size + 1;
                    row_end = i*block_size;
                    col_start = (j-1)*block_size + 1;
                    col_end = j*block_size;
                    block = saved_omega_spectral_highres(row_start:row_end, col_start:col_end, k);
                    saved_omega_highres_downsampled(i,j,k) = mean(block(:), 'omitnan');
                end
            end
        end
        
        % Calculate RMSE
        rmse_fvm = zeros(1, length(save_indices));
        for k = 1:length(save_indices)
            diff = saved_omega(:,:,k) - saved_omega_highres_downsampled(:,:,k);
            rmse_fvm(k) = sqrt(mean(diff(:).^2));
        end
        
        % Save RMSE for IC1
        save(fullfile('data', sprintf('fvm_ic1_%dx%d.mat', m, m)), 'rmse_fvm');
    end
end

% Save aggregate timing results
save(fullfile('data', 'fvm_timing_results.mat'), 'fvm_timing_results');
fprintf('Total simulation time: %.2f seconds\n', toc(total_tic));