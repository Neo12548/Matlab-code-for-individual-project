clear all;

%Set the range of frequency
fs = 0; fe = 2; 
%Set the FWHM - full width at half maximum
FWHM_f = 300e9;

Np = 1e5;
c = 3e8;
N = 203;
%lamda = 14.9e-6;
lamda = 14.6e-6;
L = N * lamda;

%Set the modal gain g0
g0 = 2.25 / L;
alpha_h = 4;


%Centre frequency in THz
fB = 2.65;
fC = 2.9;

frequency_B = fs*fB:(fe-fs)*fB/Np:fe*fB;
frequency_C = fs*fC:(fe-fs)*fC/Np:fe*fC;

%Set the refractive indices
n_offset = 3.68;

%Calculate the gain profile and then nKK
for h = 1:length(frequency_B)
fn = frequency_B(h);
g_f_B(h) = g0 * (FWHM_f/2)^2 / ((fn * 1e12 - fB * 1e12)^2 + (FWHM_f/2)^2);
n_im_B(h) = (-1/2) * g_f_B(h) * c / (2 * pi * fn * 1e12);
delta_n_re_B(h) = - alpha_h * (fn - fB) * 1e12 / FWHM_f * n_im_B(h);
n_kk_B(h) = n_offset + delta_n_re_B(h);
end

for h = 1:length(frequency_C)
fn = frequency_C(h);
g_f_C(h) = g0 * (FWHM_f/2)^2 / ((fn * 1e12 - fC * 1e12)^2 + (FWHM_f/2)^2);
n_im_C(h) = (-1/2) * g_f_C(h) * c / (2 * pi * fn * 1e12);
delta_n_re_C(h) = - alpha_h * (fn - fC) * 1e12 / FWHM_f * n_im_C(h);
n_kk_C(h) = n_offset + delta_n_re_C(h);
end

%Without gain dispersion:

%Sets the refractive indices
n1 = 3.676; n2 = 3.535;   

%Calculate the interface constants S
S12=(1/(2*n1))*[n1+n2 n1-n2;n1-n2 n1+n2];
S21=(1/(2*n2))*[n2+n1 n2-n1;n2-n1 n2+n1];


for h = 1:length(frequency_B)   % Avoid the situation when                                           

fn_B = 1e12 * frequency_B(h) ;  
fn_C = 1e12 * frequency_C(h) ;

%Calculate the optical thickness
beta_FP_B_1 = 2 * pi * (fn_B / c) * n_kk_B(h) - 1i * g_f_B(h) / 2;
beta_FP_B_2 = 2 * pi * (fn_B / c) * n_kk_B(h) - 1i * g_f_B(h) / 2;

theta_B_1 = beta_FP_B_1 * (lamda / 2) ;    
theta_B_2 = beta_FP_B_2 * (lamda / 2) ;    
theta_B_2_2 = beta_FP_B_2 * (lamda) ;      
theta_B_2_3 = beta_FP_B_2 * (2 * lamda) ;     

beta_FP_C_1 = 2 * pi * (fn_C / c) * n_kk_C(h) - 1i * g_f_C(h) / 2;
beta_FP_C_2 = 2 * pi * (fn_C / c) * n_kk_C(h) - 1i * g_f_C(h) / 2;

theta_C_1 = beta_FP_C_1 * (lamda / 2) ;   
theta_C_2 = beta_FP_C_2 * (lamda / 2) ;    
theta_C_2_2 = beta_FP_C_2 * (lamda) ;    
theta_C_2_3 = beta_FP_C_2 * (2 * lamda) ;     

%Calculate the medium constant P - B
Pn1_B = [exp(-1i*theta_B_1) 0;0 exp(1i*theta_B_1)];
Pn2_B = [exp(-1i*theta_B_2) 0;0 exp(1i*theta_B_2)];
Pn3_B = [exp(-1i*theta_B_2_2) 0;0 exp(1i*theta_B_2_2)];
Pn4_B = [exp(-1i*theta_B_2_3) 0;0 exp(1i*theta_B_2_3)];

%Construct the M matrix of the whole strcture - B
M_B = (Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)*(Pn1_B*S12*Pn3_B*S21)^2*(Pn1_B*S12*Pn2_B*S21)^5 ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^4*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^3 ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^2*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^3*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^6 ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^2*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^11 ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^2*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)^2*(Pn1_B*S12*Pn2_B*S21)^5*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)^2*(Pn1_B*S12*Pn2_B*S21)^2*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^4*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^2 ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^4 ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^5*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^2 ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^3*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^8 ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)*(Pn1_B*S12*Pn3_B*S21)^2*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^12*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)*(Pn1_B*S12*Pn3_B*S21)^2*(Pn1_B*S12*Pn2_B*S21)^5 ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^7*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21) ...
    *(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^5*(Pn1_B*S12*Pn3_B*S21)*(Pn1_B*S12*Pn2_B*S21)^3 ...
    *(Pn1_B*S12*Pn4_B*S21)*Pn1_B;

%Calculate the medium constant P - C
Pn1_C = [exp(-1i*theta_C_1) 0;0 exp(1i*theta_C_1)];
Pn2_C = [exp(-1i*theta_C_2) 0;0 exp(1i*theta_C_2)];
Pn3_C = [exp(-1i*theta_C_2_2) 0;0 exp(1i*theta_C_2_2)];
Pn4_C = [exp(-1i*theta_C_2_3) 0;0 exp(1i*theta_C_2_3)];

%Construct the M matrix of the whole strcture - C
M_C = (Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)*(Pn1_C*S12*Pn3_C*S21)^2*(Pn1_C*S12*Pn2_C*S21)^5 ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^4*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^3 ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^2*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^3*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^6 ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^2*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^11 ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^2*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)^2*(Pn1_C*S12*Pn2_C*S21)^5*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)^2*(Pn1_C*S12*Pn2_C*S21)^2*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^4*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^2 ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^4 ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^5*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^2 ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^3*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^8 ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)*(Pn1_C*S12*Pn3_C*S21)^2*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^12*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)*(Pn1_C*S12*Pn3_C*S21)^2*(Pn1_C*S12*Pn2_C*S21)^5 ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^7*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21) ...
    *(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^5*(Pn1_C*S12*Pn3_C*S21)*(Pn1_C*S12*Pn2_C*S21)^3 ...
    *(Pn1_C*S12*Pn4_C*S21)*Pn1_C;

%Calculate the reflection coefficient
Ramp_B = M_B(2,1) / M_B(1,1);
Ramp_C = M_C(2,1) / M_C(1,1);


%Calculate the power reflection coefficient
Rpow_B = (abs(Ramp_B)).^2;
Reflectance_B(h)=Rpow_B;
Rpow_C = (abs(Ramp_C)).^2;
Reflectance_C(h)=Rpow_C;

end



figure;
plot(frequency_B,Reflectance_B,'r','LineWidth',1);
hold on;
plot(frequency_C,Reflectance_C,'r','LineWidth',1);
hold off;
xlim([2.3 3.2]);

xlabel('Frequency THz');                                                   
ylabel('Power reflectivity');



