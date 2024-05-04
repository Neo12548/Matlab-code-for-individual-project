clear all;
close all;

% Initialization
fs = 0.92; fe = 1.08; 
Np = 1000;
fB = 2.9077e+12;

g = linspace(0, 0.06, Np);

n1 = 3.7; n2 = 3.5907; N = 50;
n_ave = 3.605;
tao_0 = N / (2 * fB);
delta_f = (fe - fs) / 1e6;
Normalized_frequency = linspace(fs, fe, Np+1);

% Precompute constants
S12 = (1 / (2 * n1)) * [n1 + n2, n1 - n2; n1 - n2, n1 + n2];
S21 = (1 / (2 * n2)) * [n2 + n1, n2 - n1; n2 - n1, n2 + n1];


% Calculate phase differences
Z = computePhaseDiff(Normalized_frequency, delta_f, N, Np, S12, S21, g);


% Test with n2=n1
n2 = n1;
S12 = (1 / (2 * n1)) * [n1 + n2, n1 - n2; n1 - n2, n1 + n2];
S21 = (1 / (2 * n2)) * [n2 + n1, n2 - n1; n2 - n1, n2 + n1];



Z0 = computePhaseDiff(Normalized_frequency, delta_f, N, Np, S12, S21, g);


%figure;
%contourf(Normalized_frequency, g*200, (Z ./ Z0), [0.01, 0.1, 1, 10, 100, 1000, 10000], 'LineColor', 'red','ShowText','on');
%xlabel('Normalized Frequency (f/f0)');
%ylabel('Gain');
%title('Contour Plot of Transmission Phase');

figure;
contourf(Normalized_frequency, g*200, tao_0 * (Z ./ Z0) * 1e12, [1, 2.1, 4.6, 10, 21, 46, 100], 'LineColor','red','ShowText','on');
xlabel('Normalized Frequency (f/f0)');
ylabel('Gain');
title('uniform Group delay (ps)');

%figure(3);
%contourf(Normalized_frequency, g*200, n_ave * (Z ./ Z0), [0.01, 0.1, 1, 10, 100, 1000, 10000], 'LineColor', 'red','ShowText','on');
%xlabel('Normalized Frequency (f/f0)');
%ylabel('Gain');
%title('uniform Group index');

% Function to compute phase difference
function phase_diff = computePhaseDiff(Normalized_frequency, delta_f, N, Np, S12, S21, g)
    phase_diff = zeros(1, length(Normalized_frequency));
  for n = 1:Np
    for i = 1:length(Normalized_frequency)
        fn = Normalized_frequency(i);
        theta = (pi / 2) * [fn, fn + delta_f] + 1i * g(n);
        for j = 1:2
            Pn = [exp(1i * theta(j)), 0; 0, exp(-1i * theta(j))];
            M = (Pn * S12 * Pn * S21)^N * Pn;
            Rtrans = 1 / M(1, 1);
            phase_diff(j, i) = angle(Rtrans);
        end

    phase_diff(n,i) = (phase_diff(1, i) - phase_diff(2, i)) / (2 * pi * delta_f);
    end
   
  end
    
end
