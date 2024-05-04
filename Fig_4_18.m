% Constants
c = 3e8; % Speed of light in m/s

Np = 1e6;


f01 = 2.85e12; % Center frequency in Hz
f02 = 2.95e12; % Center frequency in Hz

fs = 0; fe = 2; 
frequency_1 = fs*f01:(fe-fs)*f01/Np:fe*f01;
frequency_2 = fs*f02:(fe-fs)*f02/Np:fe*f02;

Delta_f = 300e9; % FWHM in Hz

gainL=1; 
Lc=5e-3; % cavity length
gain = gainL/(Lc); % Peak gain of 1 at the center frequency

alpha_H = 4; % Henry's factor
nr = 3.867; % Refractive index offset

% Calculate the Lorentzian gain profile and nKK
[g1, n_KK1] = calculateLandR(frequency_1, gain, Delta_f/2, f01, c, alpha_H, Delta_f, nr);
[g2, n_KK2] = calculateLandR(frequency_2, gain, Delta_f/2, f02, c, alpha_H, Delta_f, nr);

n0 = 1; n1 = nr; n2 = nr; 

S01 = (1/(2*n0))*[n0+n1 n0-n1;
                  n0-n1 n0+n1];
S10 = (1/(2*n1))*[n1+n0 n1-n0;
                  n1-n0 n1+n0];
S02 = (1/(2*n0))*[n0+n2 n0-n2;
                  n0-n2 n0+n2];
S20 = (1/(2*n2))*[n2+n0 n2-n0;
                  n2-n0 n2+n0];

% Calculate the phase and Z
[P1, phase_1, Z1, M_1, Rtrans_1] = calculatePhaseAndTransmissions(frequency_1, n_KK1, g1, S01, S10, c, Lc, Np);
[P2, phase_2, Z2, M_2, Rtrans_2] = calculatePhaseAndTransmissions(frequency_2, n_KK2, g2, S02, S20, c, Lc, Np);



figure;
yyaxis left; %Plot data on the left Y axis
plot(frequency_1/1e12, P1, 'k', 'LineWidth', 1);
hold on;
plot(frequency_2/1e12, P2, 'r', 'LineWidth', 1);
hold off;
set(gca, 'YScale', 'log');
ylim([10^0 10^3]);
ylabel('Power reflectivity');

yyaxis right; %Plot data on the right Y axis
plot(frequency_1/1e12, Z1/pi, 'k', 'LineWidth', 1);
hold on;
plot(frequency_2/1e12, Z2/pi, 'r', 'LineWidth', 1); 
hold off;
ylabel('Transmission phase (pi)');


xlim([2.823 2.837]); %10th and 11th
xlabel('Frequency THz');                                                   
title('Reflectance and Phase Spectrum of Fabry-Perot Cavity'); 



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

function [P1, phase_1, Z1, M_1, Rtrans_1] = calculatePhaseAndTransmissions(frequency_1, n_KK1, g1, S01, S10, c, Lc, Np)
    % Initializing phase_1 and Z1 matrices based on the size of frequency_1 and Np
    phase_1 = zeros(1, length(frequency_1));
    Z1 = zeros(1, Np+1);
    
    % Loop over each frequency
    for i = 1:length(frequency_1)
        beta_FP1 = 2 * pi * (frequency_1(i) / c) * n_KK1(i) - 1j * g1(i) / 2;
        theta_1 = beta_FP1 * Lc; 

        Pn1 = [exp(-1j*theta_1) 0; 0 exp(1j*theta_1)];
        
        M_1 = S01*Pn1*S10;
        Rtrans_1 = 1 / M_1(1,1);
        phase_1(i) = angle(Rtrans_1);

        Ramp_1 = M_1(2,1) / M_1(1,1);
        Rpow_1 = (abs(Ramp_1))^2;
        P1(i) = Rpow_1;
    end     
        % Nested loop for calculating Z1
    for j = 2:Np
            if j == 2
                Z1(j) = abs(phase_1(j));
            else    
                diff1 = abs(abs(phase_1(j)) - abs(phase_1(j-1)));
                Z1(j) = Z1(1, j-1) + diff1;
            end
    end
   
end

