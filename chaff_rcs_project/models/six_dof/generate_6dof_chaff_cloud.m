function [positions, orientations, lengths] = generate_6dof_chaff_cloud(N, dt, maxTime, V, radius, rho, f, a1, g)

% Constants
S = 4 * pi * radius^2;
cXPrime = 0.52; cYPrime = 0.82;
[cX, cY] = deal(cXPrime * (1 - a1), cYPrime * f * (1 - a1));
m = 0.09; % mass of chaff in kg
theta = pi; PsiS = pi;
zeta = pi; psi = pi;

positions = zeros(N, 3);
orientations = zeros(N, 3);
lengths = radius + 0.01 * randn(N,1); % randomized lengths

for i = 1:N
    x = rand()*5; y = rand()*5; z = rand()*5;
    V_local = V;

    for t = 0:dt:maxTime
        n_d = [(-sin(zeta)*cos(psi)), cos(zeta), sin(zeta)*sin(psi)];
        n_v = [(cos(theta)*cos(psi)), sin(theta), (-cos(theta)*sin(psi))];
        n_y = cross(n_v, cross(n_v,n_d) * dot(n_v,n_d)); 
        n_y = n_y / norm(n_y);

        gamma_s = -asin(((n_v(3)*n_y(1)) - (n_v(1)*n_y(3))*n_y(2)) / ...
            (sqrt(n_v(1)^2 + n_v(3)^2) * abs(n_y(2))));

        X = 0.5 * cX * rho * V_local^2 * S;
        Y = 0.5 * cY * rho * V_local^2 * S;

        dV = ((-X - m*g*sin(theta))/m)*dt;
        dTheta = ((Y*cos(gamma_s) - m*g*cos(theta))/(m*V_local))*dt;
        dPsiS = ((-Y*sin(gamma_s))/(m*V_local*cos(theta)))*dt;

        V_local = max(1, min(300, V_local + dV));
        theta = theta + dTheta;
        PsiS = PsiS + dPsiS;

        dx = V_local * cos(theta) * cos(PsiS) * dt;
        dy = V_local * sin(theta) * dt;
        dz = -V_local * cos(theta) * sin(PsiS) * dt;

        x = x + dx;
        y = y + dy;
        z = z + dz;
    end

    positions(i,:) = [x, y, z];
    vec = [cos(theta)*cos(PsiS), sin(theta), -cos(theta)*sin(PsiS)];
    orientations(i,:) = vec / norm(vec);
end

% Radar polarization vector
pol_vec = [0 1 0]; 
cos_thetas = max(min(sum(orientations .* pol_vec, 2), 1), -1);
effectiveness = 1 - cos_thetas.^2;

% Filter out NaNs
valid = all(~isnan(positions), 2) & all(~isnan(orientations), 2);
positions = positions(valid, :);
orientations = orientations(valid, :);
lengths = lengths(valid);
effectiveness = effectiveness(valid);

% === Plot ===
figure;
scatter3(positions(:,1), positions(:,2), positions(:,3), ...
    25, effectiveness, 'filled'); hold on;
colormap(jet); clim([0 1]); colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('6DoF Chaff Cloud (Colored by sin^2θ)');
cb = colorbar; cb.Label.String = 'Scattering Effectiveness (sin²θ)';

% Add orientation arrows scaled by length
scale = 1;
quiver3(positions(:,1), positions(:,2), positions(:,3), ...
    orientations(:,1) .* lengths * scale, ...
    orientations(:,2) .* lengths * scale, ...
    orientations(:,3) .* lengths * scale, ...
    0, 'k');

axis equal; view(3); grid on;

% Save output
save_figure_to_figures('cloud_sixdof_visual');

end


