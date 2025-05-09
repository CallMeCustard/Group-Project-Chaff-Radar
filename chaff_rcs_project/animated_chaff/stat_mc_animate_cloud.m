% ANIMATE_CHAFF_CLOUD
% Visualize a tumbling and spreading chaff cloud in 3D and save to GIF
% Points are colored by effectiveness; arrows show dipole orientation + length

clear; clf;

%% Parameters
N = 300;
T = 50;
lambda = 0.032;
L_mean = lambda / 2;
L_std = L_mean * 0.2;
spread_0 = 0.5;
spread_rate = 0.05;
drift_strength = 0.03;
pol_vec = [0 1 0];  % Radar polarisation vector

% Generate dipole lengths
L = L_mean + L_std * randn(N, 1);
L = max(L, lambda / 4);  % Clamp min length

%% Output GIF path
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
gif_filename = fullfile('figures', ['chaff_cloud_' timestamp '.gif']);

%% Initial positions and orientations
positions = spread_0 * randn(N, 3);
orientations = repmat([0 1 0], N, 1);

%% Fixed axis limits
max_spread = spread_0 + (T - 1) * spread_rate;
buffer = 1;
lim = max_spread + buffer;

%% Set up figure
figure;
axis([-lim lim -lim lim -lim lim]);
axis vis3d;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
hold on;

%% Animate over time
for t = 1:T
    % Update spread + positions
    spread_sigma = spread_0 + (t - 1) * spread_rate;
    positions = spread_sigma * randn(N, 3);

    % Orientation drift
    noise = drift_strength * randn(N, 3);
    orientations = orientations + noise;
    orientations = orientations ./ vecnorm(orientations, 2, 2);

    % Scattering effectiveness
    cos_thetas = sum(orientations .* pol_vec, 2);
    cos_thetas = max(min(cos_thetas, 1), -1);
    effectiveness = 1 - cos_thetas.^2;

    % Scale arrows by dipole length
    arrow_vecs = orientations .* L;

    % === PLOT ===
    cla;

    % Colored chaff points (effectiveness)
    scatter3(positions(:,1), positions(:,2), positions(:,3), ...
             25, effectiveness, 'filled');

    % Red arrows = dipole orientation + length
    quiver3(positions(:,1), positions(:,2), positions(:,3), ...
            arrow_vecs(:,1), arrow_vecs(:,2), arrow_vecs(:,3), ...
            1, 'r');

    % Colorbar
    colormap(jet);
    clim([0 1]);
    cb = colorbar;
    cb.Label.String = 'Scattering Effectiveness (sin²θ)';

    title(sprintf('Chaff Cloud (Timestep %d)', t));
    drawnow;

    % === Save GIF frame ===
    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if t == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end









