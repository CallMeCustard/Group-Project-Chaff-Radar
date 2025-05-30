function [time, rcs_over_time] = stat_mc_simulate_drift(N, lambda, L_mean, L_std, ...
    spread_sigma, T, drift_strength, radar_vec, pol_vec)
% SIMULATE_CHAFF_RCS_DRIFT
% Simulates RCS over time with dipoles starting in a fixed orientation
% and drifting (tumbling) gradually at each timestep.
%
% Inputs:
%   N              - Number of dipoles
%   lambda         - Radar wavelength (m)
%   L_mean         - Mean dipole length
%   L_std          - Std dev of dipole length
%   spread_sigma   - Spatial spread (constant for now)
%   T              - Number of time steps
%   drift_strength - How much orientation changes per step (small = slow drift)
%   radar_vec      - Radar direction vector
%   pol_vec        - Radar polarisation vector
%
% Outputs:
%   time           - Time vector (1:T)
%   rcs_over_time  - Total RCS at each timestep

% Normalize radar and polarisation vectors
radar_vec = radar_vec / norm(radar_vec);
pol_vec   = pol_vec / norm(pol_vec);

% Step 1: Initial dipole orientations — aligned horizontally (e.g., y-axis)
orientations = repmat([0 1 0], N, 1);  % All dipoles start horizontal

% Step 2: Fixed dipole lengths (can randomize once)
L = L_mean + L_std * randn(N, 1);
L = max(L, lambda / 10);  % Ensure physical lengths

% Preallocate output
rcs_over_time = zeros(T, 1);
time = 1:T;

% Loop over time
for t = 1:T
    % Orientation drift: small random 3D noise added to each orientation
    noise = drift_strength * randn(N, 3);
    orientations = orientations + noise;
    orientations = orientations ./ vecnorm(orientations, 2, 2);  % Re-normalize

    % Compute angle to radar polarisation
    cos_theta = orientations * pol_vec;
    sin_theta_sq = 1 - cos_theta.^2;

    % Compute individual dipole RCS
    sigma_i = (pi^2 .* (L.^2) ./ lambda^2) .* sin_theta_sq;

    % Total RCS
    rcs_over_time(t) = sum(sigma_i);
end

end
