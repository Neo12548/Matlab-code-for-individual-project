clear all;

close all;


fs = 0.9034; fe = 1.1034; 

% fs = 0.; fe = 1.9; 
Np = 1000;

fB=2.9077e+12;


Normalized_frequency = fs:(fe-fs)/Np:fe;
delta_f = (fe-fs)/1e6;

%Set the gain coefficient g
g=0;



%Sets the refractive indices and layer numbefr
n1=3.7; n2=3.4837; N=50;

% n1=1; n2=2; N=5;

n_ave=3.605;

% n_ave=1;

tao_0=N/(2*fB);

Z = zeros(1, Np+1); 


%Calculate the interface constants S
S12=(1/(2*n1))*[n1+n2 n1-n2;
                n1-n2 n1+n2];
S21=(1/(2*n2))*[n2+n1 n2-n1;
                n2-n1 n2+n1];





for i = 1:length(Normalized_frequency) 
fn=Normalized_frequency(i); 

theta1 = (pi/2)*fn+1i*g; %For slot lamda/4
theta2 = (pi/2)*fn+1i*g; %For lamda lamda/4
%Calculate the medium constant P
Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2 = [exp(1i*theta2) 0;0 exp(-1i*theta2)];

% Uniform Grating
M = (Pn1*S12*Pn2*S21)^N*Pn1; 


Rtrans = 1 / M(1,1);

phase_1(i) = angle(Rtrans);
end




for i = 1:length(Normalized_frequency) 
fn=Normalized_frequency(i); 

theta1 = (pi/2)*(fn+delta_f)+1i*g; %For slot lamda/4
theta2 = (pi/2)*(fn+delta_f)+1i*g; %For lamda lamda/4

Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2 = [exp(1i*theta2) 0;0 exp(-1i*theta2)];

% Uniform Grating
M = (Pn1*S12*Pn2*S21)^N*Pn1;



Rtrans = 1 / M(1,1);
 
phase_2(i) = angle(Rtrans);
end



for i = 2:Np
  
    diff = abs(abs(phase_2(i)) - abs(phase_1(i-1)));
 
    Z(1, i) = Z(1, i-1) + diff;

end





figure ;
plot(Normalized_frequency,Z, 'b', 'LineWidth',2);
xlim([fs fe])
xlabel('Normalized_frequency (f/f0)'); 
ylabel('Phase angle (pi)');

title('Transmission Phase Spectrum');


