clear all;

c = 3e8; % Speed of light in vacuum (m/s)



no = 1.0;


n1 = 2.4;
n2 = n1-0.1;


ns = 1.52;



lamda0=1550e-9;
Ls=1430e-9;Le=1700e-9;
Np=10000;  % set the sample frequency
x=1:Np;
lamda=Ls+x*(Le-Ls)/Np;   % array

f = c ./ lamda;
f_normalized = f / (c / lamda0);


d1a=lamda0/(4*n1);   %1
d1b=3*lamda0/(4*n1);   %2
d1c=4*lamda0/(4*n1);   %2.5
d1d=6*lamda0/(4*n1);   %3.5


d2=lamda0/(4*n2);




theta1a=2*pi*n1*d1a./lamda;   % array
theta1b=2*pi*n1*d1b./lamda;
theta1c=2*pi*n1*d1c./lamda;
theta1d=2*pi*n1*d1d./lamda;



theta2=2*pi*n2*d2./lamda;   % array



C1a = cos(theta1a);
S1a = sin(theta1a);
C1b = cos(theta1b);
S1b = sin(theta1b);
C1c = cos(theta1c);
S1c = sin(theta1c);
C1d = cos(theta1d);
S1d = sin(theta1d);


C2 = cos(theta2);
S2 = sin(theta2);


Do = [ 1, 1; no, -no];
Do_inv = inv(Do);
D1 = [ 1, 1; n1, -n1];
D1_inv = inv(D1);
D2 = [ 1, 1; n2, -n2];
D2_inv = inv(D2);
Ds = [ 1, 1; ns, -ns];
Ds_inv = inv(Ds);


So1 = Do_inv*D1;
S12 = D1_inv*D2;
S2s = D2_inv*Ds;

S21 = D2_inv*D1;


% Create a new 1x100 matrix for final
FFF = zeros(1, Np);


for k = 1:Np

thetaP1a = [ C1a(1,k)+S1a(1,k)*j, 0; 
            0, (C1a(1,k)-S1a(1,k)*j)];


thetaP1b = [ C1b(1,k)+S1b(1,k)*j, 0; 
            0, (C1b(1,k)-S1b(1,k)*j)];

thetaP1c = [ C1c(1,k)+S1c(1,k)*j, 0; 
            0, (C1c(1,k)-S1c(1,k)*j)];

thetaP1d = [ C1d(1,k)+S1d(1,k)*j, 0; 
            0, (C1d(1,k)-S1d(1,k)*j)];



thetaP2 = [ C2(1,k)+S2(1,k)*j, 0; 
         0, (C2(1,k)-S2(1,k)*j)];



M = ((thetaP1a*S12*thetaP2*S21)^13)*(thetaP1b*S12*thetaP2*S21)*(thetaP1a*S12*thetaP2*S21)*...
(thetaP1b*S12*thetaP2*S21)*(thetaP1d*S12*thetaP2*S21)*(thetaP1b*S12*thetaP2*S21)*...
(thetaP1a*S12*thetaP2*S21)*(thetaP1b*S12*thetaP2*S21)*((thetaP1a*S12*thetaP2*S21)^3)*...
(thetaP1c*S12*thetaP2*S21)*((thetaP1b*S12*thetaP2*S21)^2)*((thetaP1a*S12*thetaP2*S21)^13);




Ramp= M(2,1)/M(1,1);
Rpow=(abs(Ramp))^2;
FFF(1, k) = Rpow;  % Assign an element to the new matrix

end


figure
plot(f_normalized,FFF*100);
xlabel('Normalized Frequency (f/f_0)');
ylabel('Reflectivity %');
title('Simulated Reflectivity Spectrum of the Photonic Structure');







