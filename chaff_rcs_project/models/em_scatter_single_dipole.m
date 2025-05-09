function [sigma, theta_deg, effectiveness] = em_scatter_single_dipole(L, orientation, lambda, pol_vec)
% EM_SCATTER_SINGLE_DIPOLE
% Computes the RCS, orientation angle, and scattering effectiveness
% for a single dipole interacting with a radar signal.
%
% Inputs:
%   L           - Dipole length (in meters)
%   orientation - 1x3 or 3x1 dipole orientation vector
%   lambda      - Radar wavelength (in meters)
%   pol_vec     - Radar polarisation vector (unit, 1x3 or 3x1)
%
% Outputs:
%   sigma         - Radar cross section (m^2)
%   theta_deg     - Angle between dipole and polarisation (degrees)
%   effectiveness - Normalized scattering score [0, 1] (sin²θ)

% Ensure column vectors
orientation = orientation(:);
pol_vec = pol_vec(:);

% Normalize vectors
orientation = orientation / norm(orientation);
pol_vec = pol_vec / norm(pol_vec);

% Compute angle between dipole and polarisation
cos_theta = dot(orientation, pol_vec);
cos_theta = max(min(cos_theta, 1), -1);  % Clamp to avoid domain error
theta_rad = acos(cos_theta);
theta_deg = rad2deg(theta_rad);

% Scattering effectiveness = sin²θ
effectiveness = 1 - cos_theta^2;

% RCS formula (resonant dipole approximation)
sigma = (pi^2 * L^2 / lambda^2) * effectiveness;
end
