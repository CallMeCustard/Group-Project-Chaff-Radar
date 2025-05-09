%% STATIC RADAR–AIRCRAFT–CHAFF CLOUD VISUAL

clear; clf;

%% Scene Setup
lambda = 0.032;
N = 300;  % Number of dipoles in chaff cloud

% Positions (you can adjust)
aircraft_pos = [0, 0, 0];            % Aircraft at origin
radar_pos    = [0, -60, 0];          % Radar behind aircraft
cloud_center = [0, -30, 0];          % Chaff released mid-path

% Generate chaff cloud around cloud_center
spread = 2.5;
positions = cloud_center + spread * randn(N, 3);
orientations = randn(N, 3);
orientations = orientations ./ vecnorm(orientations, 2, 2);

% Compute effectiveness for coloring
pol_vec = [0 1 0];  % Radar polarisation
cos_thetas = sum(orientations .* pol_vec, 2);
cos_thetas = max(min(cos_thetas, 1), -1);
effectiveness = 1 - cos_thetas.^2;

%% Plot
figure; hold on; axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Static Radar–Aircraft–Chaff Cloud Scene');
view(3);

% Plot radar and aircraft
plot3(radar_pos(1), radar_pos(2), radar_pos(3), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
text(radar_pos(1), radar_pos(2)-3, radar_pos(3), 'Radar');

plot3(aircraft_pos(1), aircraft_pos(2), aircraft_pos(3), 'b^', 'MarkerSize', 10, 'LineWidth', 2);
text(aircraft_pos(1), aircraft_pos(2)+2, aircraft_pos(3), 'Aircraft');

% Plot radar beam (line of sight)
plot3([radar_pos(1), aircraft_pos(1)], ...
      [radar_pos(2), aircraft_pos(2)], ...
      [radar_pos(3), aircraft_pos(3)], ...
      'k--', 'LineWidth', 1.5);

% Plot chaff cloud
scatter3(positions(:,1), positions(:,2), positions(:,3), ...
         25, effectiveness, 'filled');
colormap(jet); clim([0 1]); colorbar;
cb = colorbar; cb.Label.String = 'Scattering Effectiveness (sin²θ)';

% Optional: draw orientation arrows
quiver3(positions(:,1), positions(:,2), positions(:,3), ...
        orientations(:,1), orientations(:,2), orientations(:,3), ...
        0.75, 'r');

% %% --- Occlusion Analysis (Geometric LOS Interception) ---
% 
% % Radar-to-aircraft LOS vector
% beam_vec = aircraft_pos - radar_pos;
% beam_dir = beam_vec / norm(beam_vec);
% 
% % Vector from radar to each dipole
% to_dipoles = positions - radar_pos;
% 
% % Project onto LOS direction
% proj_lengths = dot(to_dipoles, repmat(beam_dir, N, 1), 2);
% proj_points = radar_pos + proj_lengths .* beam_dir;
% 
% % Perpendicular distance from dipoles to LOS
% dist_to_beam = vecnorm(positions - proj_points, 2, 2);
% 
% % Parameters
% beam_radius = 2.5;  % How wide is the radar "corridor"
% within_beam = dist_to_beam < beam_radius;
% 
% % Display occlusion score
% occlusion_score = sum(within_beam) / N;
% fprintf('Chaff cloud occlusion score: %.2f%% of dipoles within beam corridor.\n', ...
%          100 * occlusion_score);
% 
% % Re-plot dipoles to show occlusion
% scatter3(positions(within_beam,1), positions(within_beam,2), positions(within_beam,3), ...
%          30, [1 0 0], 'filled');  % Red = blocking
% scatter3(positions(~within_beam,1), positions(~within_beam,2), positions(~within_beam,3), ...
%          20, [0.4 0.4 0.4], 'filled');  % Grey = not blocking


% Save figure
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
save_path = fullfile('..', 'figures', ['static_scene_' timestamp '.png']);
saveas(gcf, save_path);
