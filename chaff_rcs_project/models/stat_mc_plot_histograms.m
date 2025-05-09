function stat_mc_plot_histograms(orientations, lengths, pol_vec)
% STAT_MC_PLOT_HISTOGRAMS
% Plot histograms for dipole orientation angle and dipole length
%
% Inputs:
%   orientations - Nx3 matrix of dipole unit vectors
%   lengths      - Nx1 vector of dipole lengths
%   pol_vec      - 1x3 or 3x1 radar polarisation vector

% Ensure normalised polarisation vector
pol_vec = pol_vec(:)' / norm(pol_vec);

% Compute angle between each orientation and polarisation vector
cos_theta = orientations * pol_vec';
cos_theta = max(min(cos_theta, 1), -1);  % Clamp to avoid acos error
angles_deg = acosd(cos_theta);  % Convert to degrees

% Plot orientation angle histogram
figure;
subplot(1,2,1);
histogram(angles_deg, 30, 'FaceColor', [0.2 0.5 0.8]);
xlabel('Angle to Polarisation (degrees)');
ylabel('Count');
title('Dipole Orientation Angles');
xlim([0 180]);
grid on;

% Plot length distribution
subplot(1,2,2);
histogram(lengths, 30, 'FaceColor', [0.8 0.4 0.2]);
xlabel('Dipole Length (m)');
ylabel('Count');
title('Dipole Length Distribution');
grid on;
end
