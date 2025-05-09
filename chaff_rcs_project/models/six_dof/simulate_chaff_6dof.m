function [positions, orientations, rcs_array] = simulate_chaff_6dof()
% SIMULATE_CHAFF_6DOF
% Kinematic chaff cloud generation based on 6DoF motion model
% Based on: Experimental and Numerical Study of Chaff Cloud Kinetic Performance (2018)

% --- Parameters ---
N = 200;                 % Number of chaff fibres
V0 = 30;                % Initial velocity magnitude (m/s)
radius = 0.03;          % Fibre radius (m)
S = 4 * pi * radius^2;  % Effective scattering area (sphere assumption)
m = 0.09;               % Mass (g)
rho = 8;                % Atmospheric density
f = 0.1;                % Proportionality constant
a1 = 1.5;               % Impact factor
cX_p = 0.52;            % Ref. coefficient
cY_p = 0.82;

% Time setup
dt = 0.01; maxT = 5; steps = round(maxT / dt);

% Orientation vectors
zeta = pi; psi = pi/2; theta = pi/2;
E_wave = [1, 0, 0];

positions = zeros(N, 3);
orientations = zeros(N, 3);
rcs_array = zeros(N, 1);

% --- Compute aerodynamic constants ---
[cX, cY] = aero_interference(a1, f, cX_p, cY_p);

% --- Main loop over fibres ---
for i = 1:N
    x = rand()*50; y = rand()*50; z = rand()*50;
    V = V0;
    theta = pi; PsiS = pi;
    for t = 0:dt:maxT
        n_d = [(-sin(zeta)*cos(psi)), cos(zeta), sin(zeta)*sin(psi)];
        n_v = [(cos(theta)*cos(psi)), sin(theta), (-cos(theta)*sin(psi))];
        n_y = cross(n_v, cross(n_v,n_d) * dot(n_v,n_d));
        n_y = n_y / norm(n_y);

        gamma_s = -asin(((n_v(3)*n_y(1)) - (n_v(1)*n_y(3))*n_y(2)) / ...
            (sqrt(n_v(1)^2 + n_v(3)^2) * abs(n_y(2))));

        X = 0.5 * cX * rho * V^2 * S;
        Y = 0.5 * cY * rho * V^2 * S;

        dV = ((-X - m*9.81*sin(theta))/m)*dt;
        dTheta = ((Y*cos(gamma_s) - m*9.81*cos(theta))/(m*V))*dt;
        dPsiS = ((-Y*sin(gamma_s))/(m*V*cos(theta)))*dt;

        V = max(1, min(300, V + dV));
        theta = theta + dTheta;
        PsiS = PsiS + dPsiS;

        dx = V * cos(theta) * cos(PsiS) * dt;
        dy = V * sin(theta) * dt;
        dz = -V * cos(theta) * sin(PsiS) * dt;

        x = x + dx; y = y + dy; z = z + dz;
    end

    pos = [x, y, z];
    positions(i, :) = pos;
    orientations(i, :) = [cos(theta)*cos(PsiS), sin(theta), -cos(theta)*sin(PsiS)];

    % RCS via geometric backscatter
    nd = pos / norm(pos);
    inc_angle = acos(dot(E_wave, nd));
    rcs_array(i) = pi * radius^2 * sin(inc_angle);
end

% --- Plot ---
figure;
scatter3(positions(:,1), positions(:,2), positions(:,3), 20, rcs_array, 'filled');
colormap(jet); colorbar;
title('6DoF Chaff Trajectory and RCS');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on; view(3);
end

function [cX, cY] = aero_interference(a1, f, cX_p, cY_p)
    cX = cX_p * (1 - a1);
    cY = cY_p * f * (1 - a1);
end
