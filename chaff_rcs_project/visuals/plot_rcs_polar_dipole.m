% PLOT_RCS_POLAR_DIPOLE
% Plots dipole RCS as a function of orientation angle to radar polarisation

clear; clf;

% Parameters
lambda = 0.032;               % Radar wavelength (m)
L = lambda / 2;               % Dipole length (resonant)
pol_vec = [0; 1; 0];          % Radar polarisation (Y-axis)

% Sweep dipole angles from 0 to 180 degrees
theta_deg = linspace(0, 180, 360);  % Full angle sweep
rcs_vals = zeros(size(theta_deg));

% Compute RCS at each angle
for i = 1:length(theta_deg)
    % Dipole orientation at angle θ in X-Y plane
    theta_rad = deg2rad(theta_deg(i));
    orientation = [cos(theta_rad); sin(theta_rad); 0];  % Unit vector in XY plane

    % Call single dipole scatter model
    [sigma, ~, ~] = em_scatter_single_dipole(L, orientation, lambda, pol_vec);
    rcs_vals(i) = sigma;
end

% Polar Plot
figure;
polarplot(deg2rad(theta_deg), rcs_vals, 'LineWidth', 2);
rlim([0, max(rcs_vals)*1.1]);
title('Dipole RCS vs. Orientation Angle');
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.RAxis.Label.String = 'RCS (m^2)';
