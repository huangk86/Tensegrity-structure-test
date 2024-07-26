clc; clear; close all;

% Load all required data
load('C_S.mat');
load('C_b.mat');
load('V_mode_Theoretical.mat');
load('response_7_28hz.mat');
 load('response_9_7hz.mat');
load('response_10_45hz.mat');

%% Specify the desired experiment choice directly in the code
%%experiment choice (1 for 7.28Hz, 2 for 9.7Hz, 3 for 10.45Hz): ');
experiment_choice =3; % Change this value to 1, 2, or 3 based on the desired experiment

% Initialize variables based on user choice
switch experiment_choice
    case 1
        A =response_7_28hz;
        fs = 90;
        lowcut = 7.0;
        highcut = 7.5;
        num_periods = 500;
        period_length = 1 / 7.28;
        newOrder2 = [4 5 6 1 2 3 7 8 9 10 11 12];
    case 2
        A = response_9_7hz;
        fs = 90;
        lowcut = 9.5;
        highcut = 9.8;
        num_periods = 300;
        period_length = 1 / 9.7;
      newOrder2 = [4 5 6 1 2 3 7 8 9 10 11 12];
    case 3
        A =response_10_45hz;
        fs = 90;
        lowcut = 10.4;
        highcut = 10.5;
        num_periods = 400;
        period_length = 1 / 10.45;
        newOrder2 = [4 5 6 1 2 3 12 7 8 9 10 11];
    otherwise
        error('Invalid choice. Please enter 1, 2, or 3.');
end

% Process the selected data
displace = A(3:end, 4:end) / 1000;

% Map motion capture nodes to connectivity matrix
nodeData = displace(1, :);
nodes = [];
for i = 1:3:length(nodeData)
    nodes = [nodes, nodeData(i:i+2)];
end
nodesMatrix = reshape(nodes, [3, 13]);
N_exper = nodesMatrix(:, newOrder2);

% Process displacement data
numNodes = 12;
numSamples = size(displace, 1);
numCoords = 3;

allNodes = zeros(numCoords, numNodes, numSamples);
allNodesMatrix = zeros(numSamples, numNodes * numCoords);

for row = 1:numSamples
    nodeData = displace(row, 1:36);
    nodesMatrix = reshape(nodeData, [3, numNodes]);
    N = nodesMatrix(:, newOrder2);
    allNodes(:, :, row) = N;
    allNodesMatrix(row, :) = N(:);
end

% Bandpass filter parameters
[b, a] = butter(4, [lowcut highcut] / (fs / 2), 'bandpass');
filtered_allNodesMatrix = zeros(size(allNodesMatrix));

for i = 1:size(allNodesMatrix, 2)
    filtered_allNodesMatrix(:, i) = filtfilt(b, a, allNodesMatrix(:, i));
end

t = (0:numSamples-1) / fs; % Time vector

% Plot FFT of filtered responses
% figure;
% hold on;
% colors = lines(36);
% 
% for i = 1:36
%     x = filtered_allNodesMatrix(:, i)';
%     [f, X] = HannGYOUTEKIYOU(t, x, fs);
%     plot(f, X, 'Color', colors(i, :));
% end
% 
% legend(arrayfun(@(x) ['Coord ' num2str(x)], 1:36, 'UniformOutput', false), 'Location', 'eastoutside');
% title('FFT of Filtered Responses for All Coordinates');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;
% hold off;

% Calculate maximum amplitudes and phases
new_displacement = filtered_allNodesMatrix';
total_time = num_periods * period_length;
num_samples = round(total_time * fs);

maxAmplitudes = zeros(1, size(new_displacement, 1));
maxPhases = zeros(1, size(new_displacement, 1));
start_sample = 1;

for i = 1:size(new_displacement, 1)
    segment = new_displacement(i, start_sample:(start_sample + num_samples - 1));
    InputSignal.y = segment;
    InputSignal.t = (0:num_samples-1) * (1/fs);

    rawSpectrum = WindowedFFT(InputSignal, 'rect');
    Amplitude = abs(rawSpectrum.Y);
    Phase = rawSpectrum.AngY;

    [MaxAmplitude, MaxIndex] = max(Amplitude);
    MaxPhase = Phase(MaxIndex);

    maxAmplitudes(i) = MaxAmplitude;
    maxPhases(i) = rawSpectrum.AngY(MaxIndex);
end

% Create complex vectors for phase plots
RespX_complex = maxAmplitudes .* exp(1j * maxPhases * pi / 180);

figure;
h = compass(RespX_complex);

phases_deg = angle(RespX_complex) * (180 / pi);
colors = [0 1 0; 0 1 0; 0 0 1; 0 0 1];

for idx = 1:length(h)
    if (phases_deg(idx) > -45 && phases_deg(idx) <= 45)
        color_idx = 1;
    elseif (phases_deg(idx) > 45 && phases_deg(idx) <= 135)
        color_idx = 2;
    elseif (phases_deg(idx) > 135 || phases_deg(idx) <= -135)
        color_idx = 3;
    else
        color_idx = 4;
    end
    set(h(idx), 'Color', colors(color_idx, :), 'LineWidth', 2, 'DisplayName', ['Col ' num2str(idx)]);
end

set(findall(gca, 'Type', 'Line'), 'LineWidth', 2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
savefig('phase_lags.fig');
print('phase_lags', '-dmeta', '-vector', '-r600');

% Adjust amplitudes based on phases
phi = maxAmplitudes;
for i = 1:length(maxPhases)
    if maxPhases(i) < -90 || (maxPhases(i) > 90 && maxPhases(i) <= 180)
        phi(i) = -phi(i);
    end
end

ex_phi = phi' / max(abs(phi'));

% MAC calculation
start_mode = 7;
end_mode = 12;
mac = zeros(1, end_mode - start_mode + 1);

for i = start_mode:end_mode
    theory_phi = V_mode(:, i) / max(abs(V_mode(:, i)));
    mac(i - start_mode + 1) = (abs(ex_phi' * theory_phi)^2) / ((ex_phi' * ex_phi) * (theory_phi' * theory_phi));
end


% Plot experimental mode shapes
N1 = N_exper;
phi_norm = ex_phi / max(abs(ex_phi));
qq = reshape(phi_norm, 3, []);
displacement_vectors2 = 0.1 * qq;
f1 = figure('Units', 'centimeters', 'Position', [30, 15, 14, 10]);
tenseg_plot123456(N1 + displacement_vectors2, C_b, C_s, f1);
tenseg_plot555(N1, C_b, C_s, f1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'TickLabelInterpreter', 'latex', 'LineWidth', 1);
savefig('test_mode.fig');
print('test_mode', '-dmeta', '-vector', '-r600');
