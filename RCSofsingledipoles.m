% This is the code for calculating the RCS of the chaff particles, modelling them as single dipoles.
% It uses the Monte Carlo distribution modelled by LD for Chaff Distribution to obtain angles and wavelength.

%% Parameters
% Define global parameters for all simulations
lambda = 10e-2;          % Wavelength in meters (example: 10 cm)
meanLength = lambda / 2; % Mean chaff length
stdDev = 0.01;           % Standard deviation of chaff lengths (1 cm)
minLength = 0.01;        % Minimum chaff length (1 cm)
numParticles = 1000;     % Number of chaff particles

% Orientation parameters
azimuthRange = [0, 2*pi]; % Azimuth angles (radians)
elevationRange = [0, pi]; % Elevation angles (radians)

% Gaussian position distribution parameters
mu = [0, 0, 10]; % Mean position (center at [0, 0, 10])
sigma = [5, 5, 2]; % Standard deviations in X, Y, Z directions

% Monte Carlo bounds for positions
xBounds = [-10, 10]; 
yBounds = [-10, 10];
zBounds = [0, 20];

%% Generate Chaff Lengths
% Create chaff lengths using a normal distribution and truncate at minLength
chaffLengths = meanLength + stdDev * randn(numParticles, 1);
chaffLengths(chaffLengths < minLength) = minLength;

% Visualise chaff length distribution
figure;
histogram(chaffLengths, 30, 'FaceColor', 'b', 'EdgeColor', 'k');
title('Chaff Length Distribution');
xlabel('Chaff Length (m)');
ylabel('Frequency');
grid on;

%% Monte Carlo Distribution with Orientation
% Generate Monte Carlo positions
xMC = xBounds(1) + (xBounds(2) - xBounds(1)) * rand(numParticles, 1);
yMC = yBounds(1) + (yBounds(2) - yBounds(1)) * rand(numParticles, 1);
zMC = zBounds(1) + (zBounds(2) - zBounds(1)) * rand(numParticles, 1);

% Generate random orientations for Monte Carlo
azimuthMC = azimuthRange(1) + (azimuthRange(2) - azimuthRange(1)) * rand(numParticles, 1);
elevationMC = elevationRange(1) + (elevationRange(2) - elevationRange(1)) * rand(numParticles, 1);

% Convert orientations to vector components
uMC = cos(azimuthMC) .* sin(elevationMC) .* chaffLengths;
vMC = sin(azimuthMC) .* sin(elevationMC) .* chaffLengths;
wMC = cos(elevationMC) .* chaffLengths;

% Visualise Monte Carlo distribution with orientation
figure;
scatter3(xMC, yMC, zMC, 10, 'filled', 'MarkerFaceColor', [0, 0.5, 1]);
hold on;
quiver3(xMC, yMC, zMC, uMC, vMC, wMC, 0.5, 'Color', 'r', 'LineWidth', 1);
hold off;
title('Monte Carlo Chaff Distribution with Orientation Vectors');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend({'Chaff Positions', 'Orientation Vectors'}, 'Location', 'Best');
grid on; axis equal;

%% Average RCS Calculation for a Single Dipole
% Define azimuth angles for the plot
azimuthAngles = linspace(0, 2*pi, 1000); % Azimuth angles from 0 to 2π radians

% RCS formula for a dipole, assuming sinusoidal variation
RCS = (0.3377 * lambda^2 + 0.029 * lambda^2 .* sin(1.9 * azimuthAngles)) / 3;
display RCS

% Normalize RCS by lambda^2
normalizedRCS = RCS / lambda^2;

% Plot sigma/lambda^2 vs. azimuth angle
figure;
plot(azimuthAngles * (180 / pi), normalizedRCS, 'b-', 'LineWidth', 2); % Convert azimuth to degrees
xlabel('Azimuth Angle \phi (degrees)');
ylabel('\sigma / \lambda^2');
title('Normalized RCS vs. Azimuth Angle');
grid on;

% Display average normalized RCS
averageRCS = mean(RCS);
disp(['Average RCS of the chaff particles: ', num2str(averageRCS), ' m^2']);
