% Constants
c = 3e8; % Speed of light in m/s

% Define the frequency range and parameters for THz
f = linspace(2.35, 3.2, 1000) * 1e12; % Frequency range from 2.75 to 3.05 THz

f01 = 2.65e12; % Center frequency in Hz
f02 = 2.9e12; % Center frequency in Hz

Delta_f = 300e9; % FWHM in Hz

gainL=1; 
Lc=5e-3; % cavity length
gain = gainL/(Lc); % Peak gain of 1 at the center frequency

alpha_H = 4; % Henry's factor
nr = 3.867; % Refractive index offset

% Calculate the Lorentzian gain profile
[g1, n_KK1] = calculateLandR(f, gain, Delta_f/2, f01, c, alpha_H, Delta_f, nr);
[g2, n_KK2] = calculateLandR(f, gain, Delta_f/2, f02, c, alpha_H, Delta_f, nr);




figure;
plot(f/1e12,g1*Lc+g2*Lc,'r','LineWidth',1);
xlim([2.35 3.2]);
ylim([0 1.8]);
xlabel('frequency f (THz)');                                                   
ylabel('Power reflectivity');




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
