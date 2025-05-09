function sigma = analytic_bloom_rcs(t, sigma_inf, k)
% ANALYTIC_BLOOM_RCS
% Computes the RCS at time t using an exponential bloom model.
%
% Inputs:
%   t         - Time vector (same length as simulation time)
%   sigma_inf - Asymptotic (max) RCS
%   k         - Bloom rate constant
%
% Output:
%   sigma     - RCS at each timestep

sigma = sigma_inf * (1 - exp(-k * t));
end
