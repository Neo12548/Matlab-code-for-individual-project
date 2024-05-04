clear all;

delta_f0 = (2.820 - 2.805) / 15;

% Initialize Reflectance cell array
Reflectance = cell(1, 15);

% Calculate and store Reflectance for each frequency
for i = 1:15
    Reflectance{i} = ReCal(2.820 - (i - 1) * delta_f0);
end

% Plot Reflectance
figure;
for i = 1:15
    plot_figure(Reflectance{i});
end



% Function for calculating the reflectance
function [Reflectance] = ReCal (f0)

%Set the FWHM - full width at half maximum
FWHM_f = 300e9;

% Set the range of frequency
fs = 0; fe = 2; 
% Centre frequency in THz
fc = 2.9097;

% Sampling points
Np = 1e5;
c = 3e8;
N = 203;
lamda = 14.6e-6;
Lg = N * lamda;


Lc = 5e-3;


frequency = fs*fc:(fe-fs)*fc/Np:fe*fc;

% Refractive indices initialisation
kLg = 8;
n_eff = 3.605;
delta_n = kLg / N * n_eff / 2;
n1 = n_eff + delta_n;
n2 = n_eff - delta_n;


% Gain profile parameters
g0 = 4 / Lc;
alpha_h = 4;


for h = 1:length(frequency)
fn = frequency(h);
g_f(h) = g0 * (FWHM_f/2)^2 / ((fn * 1e12 - f0 * 1e12)^2 + (FWHM_f/2)^2);
n_im(h) = (-1/2) * g_f(h) * c / (2 * pi * fn * 1e12);
delta_n_re(h) = - alpha_h * (fn - f0) * 1e12 / FWHM_f * n_im(h);
n_kk(h) = n_eff + delta_n_re(h);
end


% Calculate S between boundaries
S12 = (1 / (2 * n2)) * [n2 + n1 n2 - n1 ; n2 - n1 n2 + n1];
S21 = (1 / (2 * n1)) * [n1 + n2 n1 - n2 ; n1 - n2 n1 + n2];


% Calculate the whole matrix
for h = 1:length(frequency)                                              

fn = 1e12 * frequency(h) ;  
beta_1 = 2 * pi * (fn / c) * n_kk(h) - 1i * g_f(h) / 2;
beta_2 = 2 * pi * (fn / c) * n_kk(h) - 1i * g_f(h) / 2;

theta_1 = beta_1 * (lamda / 2) ;     
theta_2 = beta_2 * (lamda / 2) ;      
theta_2_2 = beta_2 * (lamda) ;     
theta_2_3 = beta_2 * (2 * lamda) ;    

% Calculate medium  P 
Pn1 = [exp(-1i*theta_1) 0;0 exp(1i*theta_1)];
Pn2 = [exp(-1i*theta_2) 0;0 exp(1i*theta_2)];
Pn3 = [exp(-1i*theta_2_2) 0;0 exp(1i*theta_2_2)];
Pn4 = [exp(-1i*theta_2_3) 0;0 exp(1i*theta_2_3)];

M = (Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^6 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^11 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^4 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^8 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^12*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^7*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn4*S21)*Pn1;


% Calculate the reflection coefficient
Rre = M(2,1) / M(1,1);

% Calculate the power reflection coefficient
Rpow = (abs(Rre)).^2;
Reflectance(h) = Rpow;


end

end

% Function for the plot
function [] = plot_figure (Reflectance)
Np = 1e5; %Sampling points
fs = 0; fe = 2; 
fc = 2.908; % Centre frequency in THz
frequency = fs*fc:(fe-fs)*fc/Np:fe*fc;
plot(frequency, Reflectance, 'k', 'LineWidth', 1);
hold on;
xlim([2.888 2.892]);
ylim([0 100])
xlabel('Frequency THz');                                                   
ylabel('Intensity');
end






