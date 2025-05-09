function sigma_total = stat_mc_rcs(N, lambda, L_mean, L_std, spread_sigma, radar_vec, pol_vec)
% COMPUTE_RCS_DISCRETE
% Estimate total RCS of a chaff cloud composed of N dipoles.
% 
% Inputs:
%   N            - Number of dipoles
%   lambda       - Radar wavelength (meters)
%   L_mean       - Mean dipole length (meters)
%   L_std        - Standard deviation of dipole length
%   spread_sigma - Std dev for 3D Gaussian position spread
%   radar_vec    - Radar direction vector (incident wave direction)
%   pol_vec      - Polarisation vector of radar (e.g., [0;1;0] for linear pol)

% Ensure vectors are normalized
radar_vec = radar_vec / norm(radar_vec);
pol_vec   = pol_vec / norm(pol_vec);

% STEP 1: Random positions (unused in RCS but included for completeness)
positions = spread_sigma * randn(N, 3);  % Gaussian cloud

% STEP 2: Random dipole orientation vectors
orientations = randn(N, 3);  % Sample from normal dist
orientations = orientations ./ vecnorm(orientations, 2, 2);  % Normalize to unit vectors

% STEP 3: Random dipole lengths from Gaussian
L = L_mean + L_std * randn(N, 1);
L = max(L, lambda / 10);  % Clip lengths to physical minimum

% STEP 4: Compute angle θ_i between dipole and polarisation
cos_theta = orientations * pol_vec;  % Dot product with pol vector
sin_theta_sq = 1 - cos_theta.^2;     % sin^2(θ)

% STEP 5: Compute individual RCS values
% Using simplified dipole RCS: sigma_i = (pi^2 * L^2 / lambda^2) * sin^2(theta)
sigma_i = (pi^2 .* (L.^2) ./ lambda^2) .* sin_theta_sq;

% STEP 6: Total RCS
sigma_total = sum(sigma_i);

end
