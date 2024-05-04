% Constants
c = 3e8; % Speed of light in m/s

% Define the frequency range and parameters for THz
f = linspace(2.75, 3.05, 1e6) * 1e12; % Frequency range from 2.75 to 3.05 THz

f01 = 2.85e12; % Center frequency in Hz
f02 = 2.95e12; % Center frequency in Hz

Delta_f = 300e9; % FWHM in Hz

gainL=1; 
Lc=5e-3; % cavity length
gain = gainL/(Lc); % Peak gain of 1 at the center frequency

alpha_H = 4; % Henry's factor
nr = 3.867; % Refractive index offset


% Calculate the Lorentzian gain profile and nKK
[g1, n_KK1] = calculateLandR(f, gain, Delta_f/2, f01, c, alpha_H, Delta_f, nr);
[g2, n_KK2] = calculateLandR(f, gain, Delta_f/2, f02, c, alpha_H, Delta_f, nr);

n0 = 1; n1 = nr; n2 = nr; 

S01 = (1/(2*n0))*[n0+n1 n0-n1;
                  n0-n1 n0+n1];
S10 = (1/(2*n1))*[n1+n0 n1-n0;
                  n1-n0 n1+n0];
S02 = (1/(2*n0))*[n0+n2 n0-n2;
                  n0-n2 n0+n2];
S20 = (1/(2*n2))*[n2+n0 n2-n0;
                  n2-n0 n2+n0];

for i = 1:length(f)

beta_FP1 = 2 * pi * (f(i) / c) * n_KK1(i) - 1j * g1(i) / 2;
theta_1 = beta_FP1 * Lc;    
Pn1 = [exp(-1j*theta_1) 0;0 exp(1j*theta_1)];

M_1 = S01*Pn1*S10;

Ramp_1 = M_1(2,1) / M_1(1,1);
Rpow_1 = (abs(Ramp_1))^2;
Z1(i) = Rpow_1;  % Assign an element to the new matrix

end

for i = 1:length(f)

beta_FP2 = 2 * pi * (f(i) / c) * n_KK2(i) - 1j * g2(i) / 2;
theta_2 = beta_FP2 * Lc  ;    
Pn2 = [exp(-1j*theta_2) 0;0 exp(1j*theta_2)];

M_2 = S01*Pn2*S10;

Ramp_2 = M_2(2,1) / M_2(1,1);
Rpow_2 = (abs(Ramp_2))^2;
Z2(i) = Rpow_2;  % Assign an element to the new matrix

end


figure;
plot(f/1e12, Z1, 'k', 'LineWidth', 1); % f/1e12 scale is important for demonstration!!!!!
hold on;
plot(f/1e12, Z2, 'r', 'LineWidth', 1);
hold off;
xlim([2.75 3.05]);
ylim([0 400]);
yticks(0: 50: 400);
xlabel('Frequency THz');                                                   
ylabel('Power reflectivity');
title('Reflectance Spectrum of Fabry-Perot Cavity'); 


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

