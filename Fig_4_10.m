clear all;

close all;


fs = 0.92; fe = 1.08; 

% fs = 0.; fe = 1.9; 
Np = 1000;

fB=2.9077e+12;


Normalized_frequency = fs:(fe-fs)/Np:fe;
delta_f = (fe-fs)/1e6;

%Set the gain coefficient g
g = linspace(0, 0.03, Np);



%Sets the refractive indices and layer numbefr
n1=3.7; n2=3.5907; N=50;

% n1=1; n2=2; N=5;

n_ave=3.605;

% n_ave=1;

tao_0=N/(2*fB);


%Calculate the interface constants S
S12=(1/(2*n1))*[n1+n2 n1-n2;
                n1-n2 n1+n2];
S21=(1/(2*n2))*[n2+n1 n2-n1;
                n2-n1 n2+n1];


for i = 1:Np


for j = 1:length(Normalized_frequency) 
fn=Normalized_frequency(j); 

theta1 = (pi/2)*fn+1i*g(i); %For slot lamda/4
theta2 = (pi/2)*fn+1i*g(i); %For lamda lamda/4


%Calculate the medium constant P
Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2 = [exp(1i*theta2) 0;0 exp(-1i*theta2)];

% Uniform Grating
M = (Pn1*S12*Pn2*S21)^N*Pn2; 


Ramp= M(2,1)/M(1,1);
Rpow=(abs(Ramp))^2;
Z(i,j) = Rpow;  % Assign an element to the new matrix

end
end


figure
contourf(Normalized_frequency, g*200, Z, [0.01, 0.1, 1, 10, 100, 1000, 10000], 'LineColor', 'red', 'showtext', 'on');
xlabel('Normalized Frequency (f/f0)');
ylabel('Gain (gL)');
title('Countour Plot of Reflectivity');
