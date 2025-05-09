%% TEST: Occlusion Score Over Time During Chaff Bloom
% Simulates a spreading chaff cloud between radar and aircraft
% Computes occlusion score (LOS interception) at each timestep

clear; clf;

%% Parameters
lambda = 0.032;
N = 300;            % Dipoles
T = 50;             % Timesteps
spread_0 = 1.0;
spread_rate = 0.15;
drift_strength = 0.03;
beam_radius = 2.5;

% Fixed positions
aircraft_pos = [0, 0, 0];
radar_pos    = [0, -60, 0];
cloud_center = [0, -30, 0];

% Init orientation and arrays
orientations = repmat([0 1 0], N, 1);
occlusion_scores = zeros(T, 1);

%% Set up plot
figure;
set(gcf, 'Position', [100, 100, 1000, 500]);

for t = 1:T
    % Update spread and random positions
    spread_sigma = spread_0 + (t - 1) * spread_rate;
    positions = cloud_center + spread_sigma * randn(N, 3);

    % Orientation drift
    noise = drift_strength * randn(N, 3);
    orientations = orientations + noise;
    orientations = orientations ./ vecnorm(orientations, 2, 2);

    %% --- Occlusion Analysis ---
    beam_vec = aircraft_pos - radar_pos;
    beam_dir = beam_vec / norm(beam_vec);
    to_dipoles = positions - radar_pos;
    proj_lengths = dot(to_dipoles, repmat(beam_dir, N, 1), 2);
    proj_points = radar_pos + proj_lengths .* beam_dir;
    dist_to_beam = vecnorm(positions - proj_points, 2, 2);
    within_beam = dist_to_beam < beam_radius;
    occlusion_scores(t) = sum(within_beam) / N;

    %% Plot Scene (Left)
    subplot(1,2,1); cla;
    scatter3(positions(~within_beam,1), positions(~within_beam,2), positions(~within_beam,3), ...
             15, [0.4 0.4 0.4], 'filled'); hold on;
    scatter3(positions(within_beam,1), positions(within_beam,2), positions(within_beam,3), ...
             20, [1 0 0], 'filled');
    plot3([radar_pos(1), aircraft_pos(1)], ...
          [radar_pos(2), aircraft_pos(2)], ...
          [radar_pos(3), aircraft_pos(3)], 'k--', 'LineWidth', 1.5);
    plot3(radar_pos(1), radar_pos(2), radar_pos(3), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
    plot3(aircraft_pos(1), aircraft_pos(2), aircraft_pos(3), 'b^', 'MarkerSize', 10, 'LineWidth', 2);
    title(sprintf('Chaff Cloud at t = %d', t));
    xlabel('X'); ylabel('Y'); zlabel('Z'); axis equal; view(3); grid on;

    %% Plot Occlusion Score (Right)
    subplot(1,2,2);
    plot(1:t, occlusion_scores(1:t), 'r-', 'LineWidth', 2);
    ylim([0 1]); xlim([1 T]); grid on;
    xlabel('Time Step'); ylabel('Occlusion Score');
    title('Occlusion Score Over Time');

    drawnow;
    frame = getframe(gcf);
    img = frame2im(frame);
    save_gif_frame(img, 'occlusion_scene', t);  % auto-saves to /figures

end

% Save final scene and occlusion plot
save_figure_to_figures('occlusion_vs_time');
