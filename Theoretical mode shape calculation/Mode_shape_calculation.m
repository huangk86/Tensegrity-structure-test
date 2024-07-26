clc; clear; close all;

%% Basic parameters for bars
rho_b = 1275.651;
Eb = 290e09;
% Define the dimensions and mass of the bar
length_s = 1.0; % Length in meters
bar_mass = 0.04232; % Mass in kg
outer_diameter = 0.01; % Outer diameter in meters
inner_diameter = 0.0076; % Inner diameter in meters
% Calculate the cross-sectional area of the bar
outer_radius = outer_diameter / 2;
inner_radius = inner_diameter / 2;
cross_sectional_area = pi * (outer_radius^2 - inner_radius^2); % Cross-sectional area in m^2
r_b1 = 5e-3;
A_b1 = cross_sectional_area;

%% Basic parameters for Kalve strings
Es = 0.608E9;
% Define the dimensions and mass of the string
length_cm = 500; % Length in cm
mass_g = 2.45; % Mass in grams
% Measured diameters of the string in mm
diameters_mm = [0.893, 0.894, 0.886];
% Convert length to meters, mass to kg, and diameter to meters
length_m = length_cm / 100;
mass_kg = mass_g / 1000;
diameters_m = diameters_mm / 1000;
% Calculate the average diameter
average_diameter_m = mean(diameters_m);
% Calculate the average radius
average_radius_m = average_diameter_m / 2;
% Calculate the cross-sectional area (circular section)
cross_sectional_area_m2 = pi * (average_radius_m^2);
% Calculate the volume of the string
volume_m3 = cross_sectional_area_m2 * length_m;
% Calculate the density of the string
density_kg_per_m3 = mass_kg / volume_m3;

rho_s = density_kg_per_m3;
r_s1 = average_radius_m;
A_s1 = cross_sectional_area_m2;

lumped = 0; % Use lumped matrix: 1 for yes, 0 for no

%% Geometric information
p = 6; % Number of bars
r = 0.5463762439930931364355749518424; % Circumscribed circle radius
h = 0.35; % Prism height
fai = 4 * pi / p; % Twisting angle

% Node coordinates
for i = 1:p
    N(:, i) = [r * cos(2 * pi * (i - 1) / p);
               r * sin(2 * pi * (i - 1) / p);
               -0.5 * h];
end
for i = p + 1:2 * p
    N(:, i) = [r * cos(2 * pi * (i - p - 1) / p + fai);
               r * sin(2 * pi * (i - p - 1) / p + fai);
               0.5 * h];
end
NodeMatrix = N;
BarConnectivity = [1, 6, 3, 8, 5, 10; 7, 12, 9, 2, 11, 4]';
StringConnectivity1 = [1:2 * p; 2:p, 1, p + 2:2 * p, p + 1]';
StringConnectivity3 = [1, 2, 3, 4, 5, 6; 12, 7, 8, 9, 10, 11]';
StringConnectivity = [StringConnectivity1; StringConnectivity3];
Cb_in = BarConnectivity;
Cs_in = StringConnectivity;

C_b = tenseg_ind2C(Cb_in, N);
C_s = tenseg_ind2C(Cs_in, N);
C = [C_b; C_s];

[ne, nn] = size(C);
%% Plot the basic configuration
tenseg_plot123456(N, C_b, C_s);

%% Natural boundary conditions
pinned_X = [];
pinned_Y = [];
pinned_Z = [];
[Ia, Ib, a, b] = tenseg_boundary(pinned_X, pinned_Y, pinned_Z, nn);

%% Grouping
gr = {(1:6); (7:12); (13:18); (19:24)}; % Number of elements in one group
Gp = tenseg_str_gp(gr, C);
H = N * C'; % Element's direction matrix
l = sqrt(diag(H' * H)); % Element lengths
l_gp = pinv(Gp) * l; % Element lengths in groups
Cell_H = mat2cell(H, 3, ones(1, size(H, 2))); % Transfer matrix H into a cell
A_1a = Ia' * kron(C', eye(3)) * blkdiag(Cell_H{:}); % Equilibrium matrix
A_1ag = A_1a * Gp; % Equilibrium matrix in group constraints
A_2a = A_1a * diag(l .^ -1); % Equilibrium matrix
A_2ag = A_2a * Gp; % Equilibrium matrix in group constraints
[U1, U2, V1, V2, S1] = tenseg_svd(A_1ag);

% External force in equilibrium design
w0 = zeros(numel(N), 1);
w0a = Ia' * w0;
index_gp = [2, 3]; % Number of groups with designed force

force_range = [0:0.1:1, 1.5:0.5:20]; % Force range from 0-1N with 0.1N steps and 1-20N with 0.5N steps
num_forces = length(force_range); % Number of forces in the range
eigen_frequencies = zeros(num_forces, 36); % Initialize matrix to store eigenfrequencies for each force

%% Pre-stress settings
fd = 10 * ones(2, 1); % Force in bar is given as -1000
I = eye(size(Gp, 2));
e_d = I(:, index_gp); % e_d is the matrix to select group of member with designed force
l_d = e_d' * l_gp; % Length of top center circular strings
qd = fd ./ l_d;
z = (e_d' * V2) \ (qd - e_d' * pinv(A_1ag) * w0a); % Self-stress coefficient

q1_gp = pinv(A_1ag) * w0a;
q1 = Gp * q1_gp;
q2_gp = V2 * z;
q2 = Gp * q2_gp;
q_gp = q1_gp + q2_gp; % Force density in group
q = q1 + q2; % Force density
t = diag(l) * q; % Force vector
t_gp = pinv(Gp) * t; % Force in group

%% Cross-sectional area calculations
index_b = find(t < 0); % Index of bars in compression
index_s = setdiff(1:ne, index_b); % Index of strings
A_b = A_b1 * ones(numel(index_b), 1);
r_b = r_b1 * ones(numel(index_b), 1);
A_s = A_s1 * ones(numel(index_s), 1); % Area of string
r_s = r_s1 * ones(numel(index_s), 1);
I3 = eye(ne);
Ind_b = I3(:, index_b); % Index matrix for bar
Ind_s = I3(:, index_s); % Index matrix for string
A = [Ind_b, Ind_s] * [A_b; A_s]; % Cross-sectional area
A_gp = pinv(Gp) * A;
radius = [Ind_b, Ind_s] * [r_b; r_s]; % Radius
r_gp = pinv(Gp) * radius;
E = [Ind_b, Ind_s] * [Eb * ones(numel(index_b), 1); Es * ones(numel(index_s), 1)];

% Members' force & rest length
l0 = E .* A .* l ./ (t + E .* A);

% Density vector
rho = zeros(ne, 1);
rho(index_b) = rho_b;
rho(index_s) = rho_s;

% Mass matrix
mass = rho .* A .* l0;
M = tenseg_mass_matrix(mass, C, lumped); % Generate mass matrix
d = 0; % Damping coefficient
D = d * 2 * max(sqrt(mass .* E .* A ./ l0)) * eye(3 * nn); % Critical damping

%% Mode analysis
N = reshape(N(:), 3, []); % Nodal coordinate matrix
K_T = tenseg_stiff_matx2(C, N(:), E, A, l0);
[V_mode, D1] = eig(Ia' * K_T * Ia, Ia' * M * Ia); % Calculate vibration mode
d1 = diag(D1); % Eigen value
omega = real(sqrt(d1)) / (2 * pi); % Natural frequencies

%% Plot vibration modes (uncomment to enable plotting)
num_plt = 7:12; 
for i = 1:numel(num_plt)
    f1 = figure('Units', 'centimeters', 'Position', [25, 15, 10, 8]);

    % Plot vibration modes
    tenseg_plot123456(N + 0.05 * max(l0) * reshape(Ia * V_mode(:, num_plt(i)), 3, []), C_b, C_s, f1);
    tenseg_plot555(N, C_b, C_s, f1);

    % Set figure properties
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);

    % Save figure as .fig format
    savefig(['Mode', num2str(num_plt(i)), '.fig']);

    % Save figure as high-resolution EMF format
    print(['Mode', num2str(num_plt(i))], '-dmeta', '-vector', '-r600'); % Set resolution to 600 and use painters for vector format
end
