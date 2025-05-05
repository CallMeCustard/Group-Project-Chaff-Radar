%Code to model 6DOF primarily using the paper: 
%https://www.researchgate.net/publication/327064480_Experimental_and_numerical_study_of_chaff_cloud_kinetic_performance_under_the_impact_of_high_speed_airflow
%implementation of RCS from: 'Absoprtion and scattering of light by small
%particles', Bohren and Huffman (1988)

function main()
    % Initialising parameters
    theta = pi/2; %pitch angle (radians)
    Psi_s = pi/3; %yaw? angle (radians)
    zeta = pi; % (radians)
    psi = pi/2; %(radians)
    E_wave = [1, 0, 0]; %Incident EM wave direction coming from x-direction

    % Unit vectors
    n_d = [(-sin(zeta)*cos(psi)), (cos(zeta)), (sin(zeta)*sin(psi))]; %unit vector of 
    n_v = [(cos(theta)*cos(psi)), sin(theta), (-cos(theta)*sin(psi))]; %unit vector of velocity 
    n_y = cross(n_v, cross(n_v,n_d) * dot(n_v,n_d)); %unit vector of lift 
    n_y = n_y / norm(n_y);

    gamma_s = -asin(((n_v(3)*n_y(1)) - (n_v(1)*n_y(3))*n_y(2)) / (sqrt((n_v(1))^2+(n_v(3))^2) * abs(n_y(2)))); % roll angle (radians)ans)

    % Constants
    m = 0.09; %Mass of chaff fibre (grams)
    g = 9.81; %accleration due to gravity (ms^2)
    rho = 8; %atmospheric density (assumed constant)
    aOne = 1.5; %Distance impact factor
    fValue = 0.1;  %Proportionality coefficient
    dt = 0.01; %Time step
    maxTime = 5; %Final time of the simulation
    numberOfChaff = 40; %Number of chaff fibres
    sizeOfChaff = 15; %Size of chaff on th e graph
    V = 30; %Size of chaff on th e graph
    radius = 0.03; %Radius of chaff
    S = 4*pi*radius^2; %area of chaff (m^2) %want to treat like a sphere

    steps = maxTime / dt;

    % Call the main function
    [inci_array, rcs_array, rcs_normalized] = kinematic_equations(gamma_s, aOne, fValue, 0.52, 0.82, dt, maxTime, m, g, numberOfChaff, rho, S, sizeOfChaff, V, steps, zeta, psi, E_wave, radius);

    % Plot RCS vs Angle
    figure;
    plot(rad2deg(inci_array), rcs_normalized, 'o');
    xlabel('Scattering Angle (degrees)');
    ylabel('Radar Cross Section (RCS)');
    title('RCS vs Scattering Angle');
    grid on;
end

% Subfunction 1
function [cX, cY] = calculateAerodynamicInterferenceImpactFactors(aOne, fValue, cXPrime, cYPrime)  %Function calculates cX and cY due to the constants f and a1 from paper
    cX = cXPrime * (1 - aOne);
    cY = cYPrime * fValue * (1 - aOne);
    %a1 is distance impact factor and f is proportionality coefficient
end

% Subfunction 2
function [inci_array, rcs_array, rcs_normalized] = kinematic_equations(gamma_s, aOne, fValue, cXPrime, cYPrime, dt, maxTime, m, g, numberOfChaff, rho, S, sizeOfChaff, V, steps, zeta, psi, E_wave, radius)

    [cX, cY] = calculateAerodynamicInterferenceImpactFactors(aOne, fValue, cXPrime, cYPrime);
    xValues = [];
    yValues = [];
    zValues = [];

    inci_array = zeros(1, numberOfChaff);
    rcs_array = zeros(1, numberOfChaff);

    for i = 1:numberOfChaff
        x = rand()*50; y = rand()*50; z = rand()*50;
        theta = pi; 
        PsiS = pi;
        V_local = V;
        zeta_i = 2*pi*rand();  % Random initial zeta for each chaff
        psi_i = 2*pi*rand();

        for t = 0:dt:maxTime
            n_d = [(-sin(zeta)*cos(psi)), cos(zeta), sin(zeta)*sin(psi)];
            n_v = [(cos(theta)*cos(psi)), sin(theta), (-cos(theta)*sin(psi))];
            n_y = cross(n_v, cross(n_v,n_d) * dot(n_v,n_d));
            n_y = n_y / norm(n_y);

            gamma_s = -asin(((n_v(3)*n_y(1)) - (n_v(1)*n_y(3))*n_y(2)) / ...
                (sqrt(n_v(1)^2 + n_v(3)^2) * abs(n_y(2))));

            X = 0.5 * cX * rho * V_local^2 * S; %Calculates drag
            Y = 0.5 * cY * rho * V_local^2 * S; %Calculates lift

            dV = ((-X - m*g*sin(theta))/m)*dt; %Change in velocity at each time step
            dTheta = ((Y*cos(gamma_s) - m*g*cos(theta))/(m*V_local))*dt; %Change in theta at each time step
            dPsiS = ((-Y*sin(gamma_s))/(m*V_local*cos(theta)))*dt; %Change in Psi at each time step

            V_local = max(1, min(300, V_local + dV));  % clamp velocity and update the value
            theta = theta + dTheta; %Update the theta value
            PsiS = PsiS + dPsiS; %update the psi value

            %Uses these new values of V, theta and Psi to find the change in the coordinate

            dx = V_local * cos(theta) * cos(PsiS) * dt;
            dy = V_local * sin(theta) * dt;
            dz = -V_local * cos(theta) * sin(PsiS) * dt;

            x = x + dx;
            y = y + dy;
            z = z + dz;
        end

        % Store coordinates
        xValues(end+1) = x;
        yValues(end+1) = y;
        zValues(end+1) = z;

        % Use final position to compute n_d direction
        positionVec = [x, y, z];
        if norm(positionVec) == 0
            n_d = [1, 0, 0]; % fallback to E_wave direction
        else
            n_d = positionVec / norm(positionVec);
        end

        % Compute scattering angle
        inci = acos((dot(E_wave, n_d) / norm(E_wave)*norm(n_d)));


        % RCS calculation
        rcs = pi * (radius)^2 * sin(inci);

        inci_array(i) = inci; %store the values for incident angle
        rcs_array(i) = rcs; %store the values for rcs

        % Normalize RCS 
        lambda = 5e-2; % Example wavelength (1.5cm) - adjust for desiredS radar
        rcs_normalized = rcs_array / (lambda^2);% Normalize by wavelength squared
    end

    % 3D scatter plot
    figure;
    scatter3(xValues, yValues, zValues, sizeOfChaff, 'filled');
    title('Chaff Trajectory');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on;
end
