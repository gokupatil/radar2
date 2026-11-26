clc;
clear;
close all;

Pt = 10e3;
Pmin = 1e-13;
Gt = 1000;
Gr = 1000;
c = 3e8;
f = 10e9;
lambda = c/f;
sigma = 1;

Rmax = (Pt*Gt*Gr*lambda^2*sigma / ...
    ((4*pi)^3*Pmin))^(1/4);

fprintf('Maximum Radar Range = %.2f km\n', Rmax/1000);

R = linspace(100,100000,500);
F = 1:0.2:20;
S = 0.1:0.5:100;

Pr = Pt*Gt*Gr*lambda^2 ./ ((4*pi)^3*R.^4);

subplot(3,2,1);
plot(R/1000,10*log10(Pr/1e-3),'b','LineWidth',2);
grid on;
title('Received Power vs Range');
xlabel('Range (km)');
ylabel('Power (dBm)');

Ptrans = Pmin*(4*pi)^3*R.^4/(Gt*Gr*lambda^2*sigma);

subplot(3,2,2);
plot(R/1000,Ptrans/1000,'r','LineWidth',2);
grid on;
title('Transmit Power vs Range');
xlabel('Range (km)');
ylabel('Power (kW)');

Gain = Pmin*(4*pi)^3*R.^4/(Pt*Gr*lambda^2*sigma);

subplot(3,2,3);
plot(R/1000,10*log10(Gain),'m','LineWidth',2);
grid on;
title('Antenna Gain vs Range');
xlabel('Range (km)');
ylabel('Gain (dB)');

lambda_f = c./(F*1e9);

Range_f = (Pt*Gt*Gr.*lambda_f.^2*sigma / ...
    ((4*pi)^3*Pmin)).^(1/4);

subplot(3,2,4);
plot(F,Range_f/1000,'g','LineWidth',2);
grid on;
title('Maximum Range vs Frequency');
xlabel('Frequency (GHz)');
ylabel('Range (km)');

Range_s = (Pt*Gt*Gr*lambda^2.*S / ...
    ((4*pi)^3*Pmin)).^(1/4);

subplot(3,2,5);
plot(S,Range_s/1000,'k','LineWidth',2);
grid on;
title('Maximum Range vs RCS');
xlabel('RCS (m^2)');
ylabel('Range (km)');

sgtitle('Radar Range Equation Analysis');
%G=(4piAe)/lam^2
%Ae=nA
