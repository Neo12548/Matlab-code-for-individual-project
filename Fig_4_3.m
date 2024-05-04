clear all;
N = 8;     % set the number N here
no = 1.0;
n1 = 2.3;
n2 = 1.38;
nL1 = n2;
nL2 = n1;
nR1 = n1;
nR2 = n2;
ns = 1.52;
lamda1=450;
% lamda2=lamda1;
lamda2=620;

lamda0=535;
nm = no;           %set nm
dm=lamda0/(4*nm);
thetaM=pi/2;       %set thetaM



theta0=5.15*pi;

Ls=300;Le=800;
Np=1000;  % set the sample frequency
x=1:Np;
lamda=Ls+x*(Le-Ls)/Np;   % array
dL1=lamda1/(4*nL1);
dL2=lamda1/(4*nL2);

dR1=lamda2/(4*nR1);
dR2=lamda2/(4*nR2);


thetaL1=2*pi*nL1*dL1./lamda;   % array
thetaL2=2*pi*nL2*dL2./lamda;   % array

thetaR1=2*pi*nR1*dR1./lamda;   % array
thetaR2=2*pi*nR2*dR2./lamda;   % array


CL1 = cos(thetaL1);
SL1 = sin(thetaL1);
CL2 = cos(thetaL2);
SL2 = sin(thetaL2);

CR1 = cos(thetaR1);
SR1 = sin(thetaR1);
CR2 = cos(thetaR2);
SR2 = sin(thetaR2);

Cm = cos(thetaM);
Sm = sin(thetaM);


Da = [ 1, 1; no, -no];
Da_inv = inv(Da);
Db = [ 1, 1; n1, -n1];
Db_inv = inv(Db);
Dc = [ 1, 1; n2, -n2];
Dc_inv = inv(Dc);
Dd = [ 1, 1; ns, -ns];
Dd_inv = inv(Dd);

Dm = [ 1, 1; nm, -nm];
Dm_inv = inv(Dm);



Sab = Da_inv*Db;
Sbc = Db_inv*Dc;
Scd = Dc_inv*Dd;

Scb = Dc_inv*Db;

Sac = Da_inv*Dc;
Sbm = Db_inv*Dm;
Smb = Dm_inv*Db;


% Create a new 1x100 matrix for final
FFF = zeros(1, Np);


for k = 1:Np

thetaPL1 = [ CL1(1,k)+SL1(1,k)*j, 0; 
         0, (CL1(1,k)-SL1(1,k)*j)];
thetaPL2 = [ CL2(1,k)+SL2(1,k)*j, 0; 
         0, (CL2(1,k)-SL2(1,k)*j)];


thetaPR1 = [ CR1(1,k)+SR1(1,k)*j, 0; 
         0, (CR1(1,k)-SR1(1,k)*j)];
thetaPR2 = [ CR2(1,k)+SR2(1,k)*j, 0; 
         0, (CR2(1,k)-SR2(1,k)*j)];



thetaPm = [ Cm+Sm*j, 0; 
         0, (Cm-Sm*j)];


M = Sac*((thetaPL1*Scb*thetaPL2*Sbc)^N)*Sbm*thetaPm*Smb*((thetaPR1*Sbc*thetaPR2*Scb)^N)*Scd;



Ramp= M(2,1)/M(1,1);
Rpow=(abs(Ramp))^2;
FFF(1, k) = Rpow;  % Assign an element to the new matrix

end


figure
plot(lamda,FFF*100);
xlabel('wavelength (nm)');
ylabel('reflectance %');