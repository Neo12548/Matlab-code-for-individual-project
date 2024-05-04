% Constants
c = 3e8; % Speed of light in m/s

% Define the frequency range and parameters for THz
f = linspace(2.75, 3.05, 1000) * 1e12; % Frequency range from 2.75 to 3.05 THz

f01 = 2.85e12; % Center frequency in Hz
f02 = 2.95e12; % Center frequency in Hz

Delta_f = 300e9; % FWHM in Hz

gainL=1; 
Lc=5e-3; % cavity length
gain = gainL/(Lc); % Peak gain of 1 at the center frequency

alpha_H = 4; % Henry's factor
nr = 3.867; % Refractive index offset

% Calculate the Lorentzian gain profile
[g1, n_KK1] = calculateLandR(f, gain, Delta_f/2, f01, c, alpha_H, Delta_f, nr);
[g2, n_KK2] = calculateLandR(f, gain, Delta_f/2, f02, c, alpha_H, Delta_f, nr);



% Plot the Lorentzian gain profile and refractive index
figure (1); % Create a new figure window
[AX, H1, H2] = plotyy(f/1e12, g1*Lc, f/1e12, g2*Lc); % Using plotyy to plot two y-axes
xlabel('Frequency (THz)');
ylabel(AX(1), 'Gain (g_0 L_c)');
ylabel(AX(2), 'Gain (g_0 L_c)');
title('Lorentzian Gain Profile in two frequencies');
grid on; % Add a grid for easier visualization

set(AX(1), 'YColor', 'black');
set(AX(1), 'YLim', [0, 1.1]); % Set the y-axis limits for the first axis
set(AX(1), 'YTick', 0:0.1:1.1); % Set custom y-axis tick marks

set(AX(2), 'YColor', 'red');
set(AX(2), 'YLim', [0, 1.1]); % Set the y-axis limits for the second axis
set(AX(2), 'YTick', 0:0.1:1.1); % Set custom y-axis tick marks

set(H1, 'LineStyle', '-', 'Color', 'black');
set(H2, 'LineStyle', '--', 'Color', 'red');

figure (2); % Create a new figure window
[AY, H1, H2] = plotyy(f/1e12, n_KK1, f/1e12, n_KK2); % Using plotyy to plot two y-axes
xlabel('Frequency (THz)');
ylabel(AY(1), 'Refractive index (n_{KK})');
ylabel(AY(2), 'Refractive index (n_{KK})');
title('Refractive Index Change in two frequencies');
grid on; % Add a grid for easier visualization

set(AY(1), 'YColor', 'black');
set(AY(1), 'YLim', [3.865, 3.869]); % Set the y-axis limits for the first axis
set(AY(1), 'YTick', 3.865:0.001:3.869); % Set custom y-axis tick marks

set(AY(2), 'YColor', 'red');
set(AY(2), 'YLim', [3.865, 3.869]); % Set the y-axis limits for the second axis
set(AY(2), 'YTick', 3.865:0.001:3.869); % Set custom y-axis tick marks

set(H1, 'LineStyle', '-', 'Color', 'black');
set(H2, 'LineStyle', '--', 'Color', 'red');


function [g, n_KK] = calculateLandR(f, A, gamma, f0, c, alpha_H, Delta_f, nr)
    % Calculate the Lorentzian line shape function g
    g = A .* ((gamma^2) ./ ((f - f0).^2 + gamma^2));

    % Calculate the imaginary part of the refractive index change
    ni_f = -1/2 * (c * g) ./ (2 * pi * f);

    % Calculate the variation of the real component of the refractive index change
    Delta_nr_f = -alpha_H * (f - f0) / Delta_f .* ni_f;

    % Calculate the real part of the refractive index
    n_KK = nr + Delta_nr_f;
end
