function E_total = em_scatter_field_sum(positions, orientations, lambda, pol_vec, r_hat)
% EM_SCATTER_FIELD_SUM
% Computes the total far-field scattered electric field at a point (e.g. radar)
% from a cloud of dipoles using a simplified amplitude + phase model.
%
% Inputs:
%   positions    - Nx3 dipole positions
%   orientations - Nx3 unit vectors for dipole orientation
%   lambda       - Radar wavelength (m)
%   pol_vec      - Radar polarisation unit vector (1x3 or 3x1)
%   r_hat        - Unit vector toward radar (e.g. [-1 0 0])
%
% Output:
%   E_total - Complex scalar representing total electric field amplitude

% Ensure vectors are unit format
pol_vec = pol_vec(:)';
r_hat = r_hat(:)';

% Constants
k = 2 * pi / lambda;

% Scattering amplitude ∝ sin(theta) between dipole and E-field
cos_thetas = sum(orientations .* pol_vec, 2);
cos_thetas = max(min(cos_thetas, 1), -1);
A = sqrt(1 - cos_thetas.^2);  % amplitude ∝ sin(theta)

% Phase delay due to dipole positions
phases = exp(1j * k * (positions * r_hat'));

% Total field (complex scalar sum)
E_total = sum(A .* phases);
end
