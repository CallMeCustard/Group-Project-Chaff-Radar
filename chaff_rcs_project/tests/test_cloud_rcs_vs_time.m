%% RCS Evolution Over Time During Chaff Bloom (EM Scattering Models)

% Parameters
lambda = 0.032;
N = 300;
T = 50;
spread_0 = 0.5;
spread_rate = 0.05;
drift_strength = 0.03;
pol_vec = [0 1 0];
r_hat = [0 1 0];  % Radar direction (same as polarisation)

% Preallocate RCS arrays
rcs_coherent = zeros(T, 1);
rcs_incoherent = zeros(T, 1);
rcs_statistical = zeros(T, 1);  % Optional

% Initial dipole orientations
orientations = repmat([0 1 0], N, 1);

for t = 1:T
    % === Spread positions over time ===
    spread_sigma = spread_0 + (t - 1) * spread_rate;
    positions = spread_sigma * randn(N, 3);

    % === Update orientations with drift ===
    noise = drift_strength * randn(N, 3);
    orientations = orientations + noise;
    orientations = orientations ./ vecnorm(orientations, 2, 2);

    % === Coherent RCS ===
    E_total = em_scatter_field_sum(positions, orientations, lambda, pol_vec, r_hat);
    rcs_coherent(t) = abs(E_total)^2;

    % === Incoherent RCS ===
    cos_thetas = sum(orientations .* pol_vec, 2);
    cos_thetas = max(min(cos_thetas, 1), -1);
    A = sqrt(1 - cos_thetas.^2);
    rcs_incoherent(t) = sum(A.^2);

    % === Statistical RCS (optional sin² model) ===
    rcs_statistical(t) = sum(sin(acos(cos_thetas)).^2);
end

% Normalize (optional for visual comparison)
rcs_coherent = rcs_coherent / max(rcs_coherent);
rcs_incoherent = rcs_incoherent / max(rcs_incoherent);
rcs_statistical = rcs_statistical / max(rcs_statistical);

%% Plot RCS over time
figure;
plot(1:T, rcs_coherent, 'b-', 'LineWidth', 1.5); hold on;
plot(1:T, rcs_incoherent, 'r--', 'LineWidth', 1.8);
plot(1:T, rcs_statistical, 'k:', 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Normalized RCS');
legend('Coherent', 'Incoherent (Power Sum)', 'Statistical sin²(θ)', 'Location', 'southeast');
title('RCS Evolution During Chaff Cloud Bloom');
grid on;
