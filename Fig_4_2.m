clear all;                      % Clears all variables from the workspace
N = 6;                          % Number of bilayers in the filter stack
no = 1.0; n1 = 2.3; n2 = 1.38; ns = 1.52;   % Refractive indices for air, first material, second material, and substrate
lamda0 = 500;                   % Central wavelength for the filter design in nm
Ls = 300; Le = 800;             % Wavelength range from 300 nm to 800 nm
Np = 1000;                      % Number of points to sample within the wavelength range
x = 1:Np; lamda = Ls + x*(Le-Ls)/Np;  % Wavelength array over the specified range
d1 = lamda0/(4*n1); d2 = lamda0/(4*n2);  % Thickness of each layer at the central wavelength
theta1 = 2*pi*n1*d1./lamda; theta2 = 2*pi*n2*d2./lamda;  % Phase shift for each layer at each wavelength

C1 = cos(theta1); S1 = sin(theta1);     % Cosine and sine of theta1
C2 = cos(theta2); S2 = sin(theta2);     % Cosine and sine of theta2

Da = [1, 1; no, -no]; Da_inv = inv(Da);  % Matrix and its inverse for air/material interface
Db = [1, 1; n1, -n1]; Db_inv = inv(Db);  % Matrix and its inverse for material1/material2 interface
Dc = [1, 1; n2, -n2]; Dc_inv = inv(Dc);  % Matrix and its inverse for material2/substrate interface
Dd = [1, 1; ns, -ns]; Dd_inv = inv(Dd);  % Matrix and its inverse for substrate/air interface

Sab = Da_inv*Db; Sbc = Db_inv*Dc; Scd = Dc_inv*Dd; Scb = Dc_inv*Db;  % Interface matrices for each interface

Z = zeros(1, Np);    % Initialize an array to hold reflectance values

% Loop over each wavelength point to calculate reflectance
for k = 1:Np
    thetaP1 = [C1(1,k)+S1(1,k)*j, 0; 
               0, (C1(1,k)-S1(1,k)*j)];  % Matrix representing phase shifts for first material
    thetaP2 = [C2(1,k)+S2(1,k)*j, 0; 
               0, (C2(1,k)-S2(1,k)*j)];  % Matrix representing phase shifts for second material

    % Calculate overall matrix for one period of the stack
    M = Sab*((thetaP1*Sbc*thetaP2*Scb)^N)*Scd;

    Ramp = M(2,1)/M(1,1);               % Calculate amplitude reflection coefficient
    Rpow = (abs(Ramp))^2;               % Calculate power reflection coefficient
    Z(1, k) = Rpow;                     % Assign calculated reflectance to array
end

% Plot the reflectance as a function of wavelength
figure
plot(lamda, Z*100);
xlabel('wavelength (nm)');
ylabel('reflectance %');
