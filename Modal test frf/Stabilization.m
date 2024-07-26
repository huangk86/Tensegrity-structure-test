
% Load data
close all;
clc;

% Open the file
fid = fopen('Modal_test_frf.txt', 'r');

% Skip the first 60 lines of metadata
for k = 1:60
    fgetl(fid);
end

% Read the numerical data starting from line 61
data = fscanf(fid, '%f %f %f', [3, Inf]);
data = data';

% Close the file
fclose(fid);

% Extract frequency, real part, and imaginary part
frequency = data(:, 1);
real_part = data(:, 2);
imag_part = data(:, 3);
complex_response = real_part + 1i * imag_part;



% Define sampling frequency
fs = 40; % Modify the sampling frequency according to your data

% Reshape complex frequency response data to a 3D array of P-by-1-by-1
FRF = reshape(complex_response, length(complex_response), 1, 1);

% Create figure and set its size
figure;
set(gcf, 'Units', 'centimeters', 'Position', [30 15 14 10]);

% Generate stabilization diagram using 'lsrf' method
modalsd(FRF, frequency, fs, 'MaxModes', 22, 'FitMethod', 'lsrf');

% Remove automatically generated legend and title
legend('off');
title('');

% Set axes properties, labels, and line width
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 2);
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Model Order', 'FontName', 'Times New Roman', 'FontSize', 12);

% Set the line width for all lines
set(findall(gca, 'Type', 'Line'), 'LineWidth', 2);

% Set right Y-axis properties (assuming no data is displayed on the right Y-axis)
yyaxis right;
ylabel('Magnitude (g/N)', 'FontName', 'Times New Roman', 'FontSize', 12); % No data plotted
xlim([6.5 12]); % Set X-axis display range from 6 Hz to 12 Hz

% Export high-resolution image
set(gca, 'FontName', 'Times New Roman');
savefig('Stabilization_Diagram.fig');
print('Stabilization_Diagram', '-dmeta', '-vector', '-r600'); % High-resolution EMF image

% Physical frequency points
phfr = [7.27, 9.7, 10.44, 11.3886];

% Perform modal fitting using 'lsrf' method
[fn, dr, ms] = modalfit(FRF, frequency, fs, 15, 'PhysFreq', phfr, 'FitMethod', 'lsrf');
