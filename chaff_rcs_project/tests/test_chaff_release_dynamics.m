%% TEST: Dynamic Chaff Release from Moving Aircraft
% Aircraft moves, chaff released over time, drifts & spreads
% Radar beam fixed; occlusion score tracked and plotted

clear; clf;

%% Parameters
lambda = 0.032;
T = 60;                          % Timesteps
dipoles_per_step = 20;
g = [0, 0, -9.8e-3];             % Simple gravity (adjust scale)
drift_strength = 0.03;
spread_std = 0.2;
beam_radius = 2.5;

w1 = 0.6; w2 = 0.4;              % Weights for effectiveness score

%% Aircraft motion
aircraft_pos = [0, 0, 0];
u_aircraft = [0, 0.5, 0];        % Forward motion in Y

%% Radar setup
radar_pos = [0, -60, 0];
beam_vec = aircraft_pos - radar_pos;
beam_dir = beam_vec / norm(beam_vec);

rcs_coherent = zeros(T,1);
rcs_incoherent = zeros(T,1);
occlusion_scores = zeros(T,1);
effectiveness = zeros(T,1);
pol_vec = [0 1 0];              % Radar polarisation

%% Arrays to hold chaff data
positions = []; velocities = []; orientations = [];  % Nx3

%% Set up figure
figure;
set(gcf, 'Position', [100, 100, 1000, 500]);

for t = 1:T
    %% Update aircraft position
    aircraft_pos = aircraft_pos + u_aircraft;

    %% Chaff release
    new_pos = repmat(aircraft_pos, dipoles_per_step, 1);
    eject_velocity = repmat([-0.2, 0, 0], dipoles_per_step, 1);  % backward
    new_vel = eject_velocity + 0.1 * randn(dipoles_per_step, 3); % random spread
    new_orient = randn(dipoles_per_step, 3);
    new_orient = new_orient ./ vecnorm(new_orient, 2, 2);

    %% Append new chaff
    positions = [positions; new_pos];
    velocities = [velocities; new_vel];
    orientations = [orientations; new_orient];

    N = size(positions, 1);

    %% Update all chaff: drift, gravity, motion
    drift = drift_strength * randn(N, 3);
    orientations = orientations + drift;
    orientations = orientations ./ vecnorm(orientations, 2, 2);

    velocities = velocities + repmat(g, N, 1);
    positions = positions + velocities;

    %% Occlusion analysis
    to_dipoles = positions - radar_pos;
    proj_lengths = dot(to_dipoles, repmat(beam_dir, N, 1), 2);
    proj_points = radar_pos + proj_lengths .* beam_dir;
    dist_to_beam = vecnorm(positions - proj_points, 2, 2);
    within_beam = dist_to_beam < beam_radius;
    occlusion_scores(t) = sum(within_beam) / N;

    %% RCS Evaluation
    if N > 0
        % Coherent RCS (field sum)
        E_total = em_scatter_field_sum(positions, orientations, lambda, pol_vec, beam_dir);
        rcs_coherent(t) = abs(E_total)^2;

        % Incoherent RCS (power sum)
        cos_thetas = sum(orientations .* pol_vec, 2);
        cos_thetas = max(min(cos_thetas, 1), -1);
        a = sqrt(1 - cos_thetas.^2);
        rcs_incoherent(t) = sum(a.^2);
    else
        rcs_coherent(t) = 0;
        rcs_incoherent(t) = 0;
    end

    %% Effectiveness score (with normalized RCS inputs)
    rcs_coh = rcs_coherent(t);
    rcs_inc = rcs_incoherent(t);
    rcs_sum = rcs_coh + rcs_inc;
    if rcs_sum > 0
        rcs_coh = rcs_coh / rcs_sum;
        rcs_inc = rcs_inc / rcs_sum;
    end
    effectiveness(t) = effectiveness_score(occlusion_scores(t), rcs_coh, rcs_inc, w1, w2);

    %% Plot scene (left)
    subplot(1,2,1); cla;
    scatter3(positions(~within_beam,1), positions(~within_beam,2), positions(~within_beam,3), ...
        8, [0.4 0.4 0.4], 'filled'); hold on;
    scatter3(positions(within_beam,1), positions(within_beam,2), positions(within_beam,3), ...
        12, [1 0 0], 'filled');

    plot3(radar_pos(1), radar_pos(2), radar_pos(3), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
    plot3(aircraft_pos(1), aircraft_pos(2), aircraft_pos(3), 'b^', 'MarkerSize', 10, 'LineWidth', 2);
    plot3([radar_pos(1), aircraft_pos(1)], ...
          [radar_pos(2), aircraft_pos(2)], ...
          [radar_pos(3), aircraft_pos(3)], 'k--', 'LineWidth', 1.5);
    xlabel('X'); ylabel('Y'); zlabel('Z'); axis equal; view(3); grid on;
    title(sprintf('Chaff Deployment at t = %d', t));

    %% Plot occlusion and RCS (right)
    subplot(1,2,2); cla;
    yyaxis left
    h_occ = plot(1:t, occlusion_scores(1:t), 'r-', 'LineWidth', 2);
    ylabel('Occlusion Score'); ylim([0 1]);

    yyaxis right
    hold on;
    h_coh = plot(1:t, rcs_coherent(1:t), 'b--', 'LineWidth', 1.5);
    h_inc = plot(1:t, rcs_incoherent(1:t), 'k:', 'LineWidth', 1.5);
    ylabel('RCS (Normalized)'); ylim auto;

    xlabel('Time Step'); xlim([1 T]); grid on;
    title('Occlusion + RCS Over Time');
    legend([h_occ, h_coh, h_inc], {'Occlusion', 'Coherent RCS', 'Incoherent RCS'}, ...
           'Location', 'northeast', 'Box', 'off');

    %% Optional: save animation frame to GIF
    frame = getframe(gcf);
    img = frame2im(frame);
    save_gif_frame(img, 'chaff_release_scene', t);
end

%% Normalize RCS values for consistent plotting
if max(rcs_coherent) > 0
    rcs_coherent = rcs_coherent / max(rcs_coherent);
end
if max(rcs_incoherent) > 0
    rcs_incoherent = rcs_incoherent / max(rcs_incoherent);
end

%% Final plot for effectiveness
figure;
plot(1:T, effectiveness, 'g-', 'LineWidth', 2);
xlabel('Time Step'); ylabel('Effectiveness Score');
ylim([0 1]); grid on;
title('Chaff Cloud Effectiveness Over Time');
save_figure_to_figures('effectiveness_score');

%% Save final static summary plot
save_figure_to_figures('chaff_release_occlusion');

