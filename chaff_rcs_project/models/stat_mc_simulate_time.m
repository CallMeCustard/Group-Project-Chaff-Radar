function [time, rcs_over_time] = stat_mc_simulate_time(N, lambda, L_mean, L_std, ...
    spread_0, spread_rate, T, radar_vec, pol_vec)
% SIMULATE_CHAFF_RCS_TIME
% Simulates RCS of a chaff cloud over T time steps as it blooms (spreads out).
%
% Inputs:
%   N            - Number of dipoles
%   lambda       - Radar wavelength (m)
%   L_mean       - Mean dipole length
%   L_std        - Std dev of dipole length
%   spread_0     - Initial cloud standard deviation (m)
%   spread_rate  - Spread growth per timestep (m/step)
%   T            - Number of time steps
%   radar_vec    - Radar propagation direction
%   pol_vec      - Radar polarisation direction
%
% Outputs:
%   time             - Time vector [1:T]
%   rcs_over_time    - Total RCS at each time step

% Normalize input vectors
radar_vec = radar_vec / norm(radar_vec);
pol_vec   = pol_vec / norm(pol_vec);

rcs_over_time = zeros(T, 1);
for t = 1:T
    spread_sigma = spread_0 + (t-1) * spread_rate;
    rcs_over_time(t) = stat_mc_rcs(N, lambda, L_mean, L_std, ...
                                            spread_sigma, radar_vec, pol_vec);
end

time = 1:T;
end
