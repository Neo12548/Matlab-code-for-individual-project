clear all;
close all;

% Initialization
fs = 0; fe = 2; 
Np = 150000;
fB = 2.9077e+12;

g=0;

c=3e8;

lamda0=69364e-5;

n1 = 3.7; n2 = 3.4837; N = 50;
%In the case of non-grating strucuture for all three structures
n1_0=3.7; n2_0=n1_0;


n_ave = 3.605;
delta_n = n1-n2;
tao_0 = N / (2 * fB);
delta_f = (fe - fs) / 1e6;
Normalized_frequency = linspace(fs, fe, Np+1);


S12 = (1 / (2 * n1)) * [n1 + n2, n1 - n2; n1 - n2, n1 + n2];
S21 = (1 / (2 * n2)) * [n2 + n1, n2 - n1; n2 - n1, n2 + n1];

S12_0=(1 / (2 * n1_0)) * [n1_0 + n2_0, n1_0 - n2_0; n1_0 - n2_0, n1_0 + n2_0];
S21_0=(1 / (2 * n2_0)) * [n2_0 + n1_0, n2_0 - n1_0; n2_0 - n1_0, n2_0+ n1_0];


for i = 1:length(Normalized_frequency) 
fn=Normalized_frequency(i); 

theta1 = (pi/2)*fn+1i*g; %For slot lamda/4
theta2 = (pi/2)*fn+1i*g; %For slot lamda/4
theta1b=3*(pi/2)*fn+1i*g; %2
theta1c=4*(pi/2)*fn+1i*g; %2.5
theta1d=6*(pi/2)*fn+1i*g; %3.5

theta2d = (pi)*fn+1i*g;  %For slot lamda/2

Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2 = [exp(1i*theta2) 0;0 exp(-1i*theta2)];
Pn1b = [exp(1i*theta1b) 0;0 exp(-1i*theta1b)];
Pn1c = [exp(1i*theta1c) 0;0 exp(-1i*theta1c)];
Pn1d = [exp(1i*theta1d) 0;0 exp(-1i*theta1d)];

Pn2d = [exp(1i*theta1d) 0;0 exp(-1i*theta1d)];

%periodic grating strcture
M1 = (Pn1*S12*Pn2*S21)^N*Pn1;
%non-grating structure for periodic
M1_0 = (Pn1*S12_0*Pn2*S21_0)^N*Pn1;  % value N can be modified at beggining

%Single-defect grating strcture
M2 = (Pn1*S12*Pn2*S21)^24*(Pn1*S12*Pn2d*S21)*(Pn1*S12*Pn2*S21)^24*Pn2; % It's a 24/1/24 structure
%non-grating structure for Single-defect
M2_0 = (Pn1*S12_0*Pn2*S21_0)^24*(Pn1*S12_0*Pn2d*S21_0)*(Pn1*S12_0*Pn2*S21_0)^24*Pn2; 


%Aperiodic grating strcture
M3 = ((Pn1*S12*Pn2*S21)^13)*(Pn1b*S12*Pn2*S21)*(Pn1*S12*Pn2*S21)*(Pn1b*S12*Pn2*S21)*...
(Pn1d*S12*Pn2*S21)*(Pn1b*S12*Pn2*S21)*(Pn1*S12*Pn2*S21)*(Pn1b*S12*Pn2*S21)*((Pn1*S12*Pn2*S21)^3)*...
(Pn1c*S12*Pn2*S21)*((Pn1b*S12*Pn2*S21)^2)*((Pn1*S12*Pn2*S21)^13)*Pn2;
%non-grating structure for Aperiodic grating strcture
M3_0 = ((Pn1*S12_0*Pn2*S21_0)^13)*(Pn1b*S12_0*Pn2*S21_0)*(Pn1*S12_0*Pn2*S21_0)*...
(Pn1b*S12_0*Pn2*S21_0)*(Pn1d*S12_0*Pn2*S21_0)*(Pn1b*S12_0*Pn2*S21_0)*(Pn1*S12_0*Pn2*S21_0)*(Pn1b*S12_0*Pn2*S21_0)*...
((Pn1*S12_0*Pn2*S21_0)^3)*(Pn1c*S12_0*Pn2*S21_0)*((Pn1b*S12_0*Pn2*S21_0)^2)*((Pn1*S12_0*Pn2*S21_0)^13)*Pn2;


%Calculate the transmission coefficient for three structures
Rtrans1 = 1 / M1(1,1);
Rtrans1_0 = 1 / M1_0(1,1);

Rtrans2 = 1 / M2(1,1);
Rtrans2_0 = 1 / M2_0(1,1);

Rtrans3 = 1 / M3(1,1);
Rtrans3_0 = 1 / M3_0(1,1);

%Calculate the phase for three structures
phase1(i) = angle(Rtrans1);
phase1_0(i) = angle(Rtrans1_0);

phase2(i) = angle(Rtrans2);
phase2_0(i) = angle(Rtrans2_0);


phase3(i) = angle(Rtrans3);
phase3_0(i) = angle(Rtrans3_0);
end


for k = 1:length(Normalized_frequency)
fn(k)=Normalized_frequency(k); 

if k==1
phase_accu1(k) = abs(phase1(k));
phase_accu1_0(k) = abs(phase1_0(k));

phase_accu2(k) = abs(phase2(k));
phase_accu2_0(k) = abs(phase2_0(k));

phase_accu3(k) = abs(phase3(k));
phase_accu3_0(k) = abs(phase3_0(k));

else

phase_accu1(k) = abs(abs(phase1(k)) - abs(phase1(k-1))) + phase_accu1(k-1);
phase_accu1_0(k) = abs(abs(phase1_0(k)) - abs(phase1_0(k-1))) + phase_accu1_0(k-1);

phase_accu2(k) = abs(abs(phase2(k)) - abs(phase2(k-1))) + phase_accu2(k-1);
phase_accu2_0(k) = abs(abs(phase2_0(k)) - abs(phase2_0(k-1))) + phase_accu2_0(k-1);

phase_accu3(k) = abs(abs(phase3(k)) - abs(phase3(k-1))) + phase_accu3(k-1);
phase_accu3_0(k) = abs(abs(phase3_0(k)) - abs(phase3_0(k-1))) + phase_accu3_0(k-1);

end 

n_eff1(k) = (c/lamda0)*(phase_accu1(k)/(2*pi*fn(k)));
n_eff1_0(k) = (c/lamda0)*(phase_accu1_0(k)/(2*pi*fn(k)));

n_eff2(k) = (c/lamda0)*(phase_accu2(k)/(2*pi*fn(k)));
n_eff2_0(k) = (c/lamda0)*(phase_accu2_0(k)/(2*pi*fn(k)));

n_eff3(k) = (c/lamda0)*(phase_accu3(k)/(2*pi*fn(k)));
n_eff3_0(k) = (c/lamda0)*(phase_accu3_0(k)/(2*pi*fn(k)));
end

figure;
% plot(Normalized_frequency,(n_ave*(n_eff1./n_eff1_0))); 
plot(Normalized_frequency,(n_ave*(phase_accu1./phase_accu1_0))); 

xlim([0.9 1.1]);
ylim([3.55 3.66]);
xlabel('Normalized_frequency (f/f0)'); 
ylabel('Effective index value');
title('Periodic structure');
%Code from 'figure(1)' to this line activate the plot for Periodic Structure

%figure(2);
%plot(Normalized_frequency,(n_ave*(n_eff2./n_eff2_0))); 
%xlim([0.9 1.1]);
%ylim([3.55 3.66]);
%xlabel('Normalized_frequency (f/f0)'); 
%ylabel('Effective index value');
%title('PhaseShift structure');
%Code from 'figure(2)' to this line activate the plot for PhaseShift Structure

%figure(3);
%plot(Normalized_frequency,(n_ave*(n_eff3./n_eff3_0))); 
%xlim([0.9 1.1]);
%ylim([3.55 3.66]);
%xlabel('Normalized_frequency (f/f0)'); 
%ylabel('Effective index value');
%title('Aperiodic structure');
%Code from 'figure(3)' to this line activate the plot for Aperiodic Structure