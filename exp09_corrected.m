clc;
clear;
close all;

T = 150;
G = 35;
Pc = -100;
B = 1e6;
k = 1.38e-23;

N = 10*log10(k*T*B);

CN = Pc - N;

GT = G - 10*log10(T);

fprintf('Noise Power = %.2f dBW\n', N);
fprintf('C/N Ratio   = %.2f dB\n', CN);
fprintf('G/T Ratio   = %.2f dB/K\n', GT);

Tr = 50:10:500;
Gr = 20:50;

Noise = 10*log10(k*Tr*B);

CN_values = Pc - Noise;

GT_gain = Gr - 10*log10(T);

GT_temp = G - 10*log10(Tr);

subplot(2,2,1);
plot(Tr,Noise,'b','LineWidth',2);
grid on;
xlabel('Temperature (K)');
ylabel('Noise Power (dBW)');
title('Noise Power vs Temperature');

subplot(2,2,2);
plot(Tr,CN_values,'r','LineWidth',2);
grid on;
xlabel('Temperature (K)');
ylabel('C/N (dB)');
title('C/N Ratio vs Temperature');

subplot(2,2,3);
plot(Gr,GT_gain,'m','LineWidth',2);
grid on;
xlabel('Gain (dBi)');
ylabel('G/T (dB/K)');
title('G/T vs Gain');

subplot(2,2,4);
plot(Tr,GT_temp,'g','LineWidth',2);
grid on;
xlabel('Temperature (K)');
ylabel('G/T (dB/K)');
title('G/T vs Temperature');

sgtitle('Link Budget and Noise Analysis');
