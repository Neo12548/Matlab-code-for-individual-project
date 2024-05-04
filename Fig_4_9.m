clear all;
close all;

% Initialization
fs = 0.92; fe = 1.08; 
Np = 30000;
fB = 2.9077e+12;

g=0;

n1 = 3.7; n2 = 3.5907; N = 50;
n_ave = 3.605;
delta_n = n1-n2;
tao_0 = N / (2 * fB);
delta_f = (fe - fs) / 1e6;
Normalized_frequency = linspace(fs, fe, Np+1);

% Precompute constants
S12 = (1 / (2 * n1)) * [n1 + n2, n1 - n2; n1 - n2, n1 + n2];
S21 = (1 / (2 * n2)) * [n2 + n1, n2 - n1; n2 - n1, n2 + n1];


% Calculate phase differences
Z1 = computePhaseDiff1(Normalized_frequency, delta_f, N, S12, S21, g);
Z2 = computePhaseDiff2(Normalized_frequency, delta_f, N, S12, S21, g);
Z3 = computePhaseDiff3(Normalized_frequency, delta_f, N, S12, S21, g);


% Test with n2=n1
n2 = n1;
S12 = (1 / (2 * n1)) * [n1 + n2, n1 - n2; n1 - n2, n1 + n2];
S21 = (1 / (2 * n2)) * [n2 + n1, n2 - n1; n2 - n1, n2 + n1];

Z01 = computePhaseDiff1(Normalized_frequency, delta_f, N, S12, S21, g);
Z02 = computePhaseDiff2(Normalized_frequency, delta_f, N, S12, S21, g);
Z03 = computePhaseDiff2(Normalized_frequency, delta_f, N, S12, S21, g);

% Plot all three on a single graph
figure(1); % Create a new figure
hold on; % Retain current plot when adding new plots

plot(Normalized_frequency, tao_0 * (Z1 ./ Z01) * 1e12, 'r', 'LineWidth', 2); % Red line for first plot
plot(Normalized_frequency, tao_0 * (Z2 ./ Z02) * 1e12, 'b', 'LineWidth', 2); % Blue line for second plot
plot(Normalized_frequency, tao_0 * (Z3 ./ Z03) * 1e12, 'g', 'LineWidth', 2); % Green line for third plot

xlim([fs fe]); % Set x-axis limits
xlabel('Normalized frequency (f/f0)'); % X-axis label
ylabel('Group delay (ps)'); % Y-axis label
legend('Periodic', 'SingleDefect', 'Aperiodic'); % Add a legend
hold off; % Release the plot to avoid further overlay

% Function to compute phase difference
function phase_diff = computePhaseDiff1(Normalized_frequency, delta_f, N, S12, S21, g)
    phase_diff = zeros(1, length(Normalized_frequency));
    for i = 1:length(Normalized_frequency)
        fn = Normalized_frequency(i);
        theta = (pi / 2) * [fn, fn + delta_f] + 1i * g;
        for j = 1:2
            Pn = [exp(1i * theta(j)), 0; 0, exp(-1i * theta(j))];
            M = (Pn * S12 * Pn * S21)^N * Pn;
            Rtrans = 1 / M(1, 1);
            phase_diff(j, i) = angle(Rtrans);
        end
    end
    phase_diff = (phase_diff(1, :) - phase_diff(2, :)) / (2 * pi * delta_f);
end

% Function to compute phase difference
function phase_diff = computePhaseDiff2(Normalized_frequency, delta_f, N, S12, S21, g)
    phase_diff = zeros(1, length(Normalized_frequency));
    for i = 1:length(Normalized_frequency)
        fn = Normalized_frequency(i);
        theta1 = (pi / 2) * [fn, fn + delta_f] + 1i * g;
        theta1d = (pi) * [fn, fn + delta_f] + 1i * g;
        for j = 1:2
            Pn1 = [exp(1i * theta1(j)), 0; 0, exp(-1i * theta1(j))];
            Pn1d = [exp(1i * theta1d(j)), 0; 0, exp(-1i * theta1d(j))];

            M = (Pn1*S12*Pn1*S21)^25*(Pn1*S12*Pn1d*S21)*(Pn1*S12*Pn1*S21)^25*Pn1; 

            Rtrans = 1 / M(1, 1);
            phase_diff(j, i) = angle(Rtrans);
        end
    end
    phase_diff = (phase_diff(1, :) - phase_diff(2, :)) / (2 * pi * delta_f);
end

% Function to compute phase difference
function phase_diff = computePhaseDiff3(Normalized_frequency, delta_f, N, S12, S21, g)
    phase_diff = zeros(1, length(Normalized_frequency));
    for i = 1:length(Normalized_frequency)
        fn = Normalized_frequency(i);
        theta1 = (pi / 2) * [fn, fn + delta_f] + 1i * g;

        theta1b=3*(pi/2)* [fn, fn + delta_f] + 1i * g; %2
        theta1c=4*(pi/2)* [fn, fn + delta_f] + 1i * g; %2.5
        theta1d=6*(pi/2)* [fn, fn + delta_f] + 1i * g; %3.5

        for j = 1:2
            Pn1 = [exp(1i * theta1(j)), 0; 0, exp(-1i * theta1(j))];

            Pn1b = [exp(1i * theta1b(j)), 0; 0, exp(-1i * theta1b(j))];
            Pn1c = [exp(1i * theta1c(j)), 0; 0, exp(-1i * theta1c(j))];
            Pn1d = [exp(1i * theta1d(j)), 0; 0, exp(-1i * theta1d(j))];


            M = ((Pn1*S12*Pn1*S21)^13)*(Pn1b*S12*Pn1*S21)*(Pn1*S12*Pn1*S21)*(Pn1b*S12*Pn1*S21)*(Pn1d*S12*Pn1*S21)*(Pn1b*S12*Pn1*S21)*(Pn1*S12*Pn1*S21)*(Pn1b*S12*Pn1*S21)*((Pn1*S12*Pn1*S21)^3)*(Pn1c*S12*Pn1*S21)*((Pn1b*S12*Pn1*S21)^2)*((Pn1*S12*Pn1*S21)^13)*Pn1;


            Rtrans = 1 / M(1, 1);
            phase_diff(j, i) = angle(Rtrans);
        end
    end
    phase_diff = (phase_diff(1, :) - phase_diff(2, :)) / (2 * pi * delta_f);
end