clear all;
fs = 0.95; fe = 1.05; 
Np = 10000;
Normalized_frequency = fs:(fe-fs)/Np:fe;

n_eff = 3.68; n1 = 3.676; n2 = 3.576;  

N = 203;

lamda = 14.9e-6;
Lc = N * lamda;

%Calculate the interface constants S
S12=(1/(2*n1))*[n1+n2 n1-n2;n1-n2 n1+n2];
S21=(1/(2*n2))*[n2+n1 n2-n1;n2-n1 n2+n1];

%Set the gain coefficient g (g_frequency = 1/2 g_power)
gL = 0;
g = gL / (4 * N);

for i = 1:length(Normalized_frequency)                                              
fn=Normalized_frequency(i);    
%Calculate the optical thickness
theta1 = (pi/2)*fn+1i*g;    
theta2 = (pi/2)*fn+1i*g;     
theta2_2 = (pi)*fn+1i*2*g;      
theta2_3 = (2*pi)*fn+1i*4*g;
%Calculate the medium constant P
Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2 = [exp(1i*theta2) 0;0 exp(-1i*theta2)];
Pn2_2 = [exp(1i*theta2_2) 0;0 exp(-1i*theta2_2)];
Pn2_3 = [exp(1i*theta2_3) 0;0 exp(-1i*theta2_3)];
%Construct the M matrix of the whole strcture
M = (Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^6 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^11 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^4 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^8 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^12*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^7*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn2_3*S21);

%Calculate the reflection coefficient
Ramp = M(2,1) / M(1,1);
%Calculate the power reflection coefficient
Rpow = (abs(Ramp)).^2;
Reflectance_1(i)=Rpow;
end



figure;
plot(Normalized_frequency,Reflectance_1,'r','LineWidth',1);
xlim([0.95 1.05]);
ylim([0 0.5]);
xlabel('Normalized_frequency (f/fB)');                                                   
ylabel('Power reflectivity');


