clc;
clear;
close all;

d = 36000;
f_u = 14;
f_d = 12;
T = 500;
k = 1.38e-23;
Latm = 2;
Lmisc = 1;
Pr_req = -100;

Pt_u = 20;
Gt_u = 45;
Gr_u = 35;
Ltx_u = 2;
Lrx_u = 1;

Pt_d = 10;
Gt_d = 35;
Gr_d = 40;
Ltx_d = 1;
Lrx_d = 2;

N0 = 10*log10(k*T);

EIRP_u = Pt_u + Gt_u - Ltx_u;
LFS_u = 92.45 + 20*log10(f_u) + 20*log10(d);
Pr_u = EIRP_u + Gr_u - LFS_u - Lrx_u - Latm - Lmisc;
CN0_u = Pr_u - N0;
LM_u = Pr_u - Pr_req;

EIRP_d = Pt_d + Gt_d - Ltx_d;
LFS_d = 92.45 + 20*log10(f_d) + 20*log10(d);
Pr_d = EIRP_d + Gr_d - LFS_d - Lrx_d - Latm - Lmisc;
CN0_d = Pr_d - N0;
LM_d = Pr_d - Pr_req;

fprintf('\n       UPLINK       DOWNLINK\n');
fprintf('EIRP   %8.2f      %8.2f\n', EIRP_u, EIRP_d);
fprintf('LFS    %8.2f      %8.2f\n', LFS_u, LFS_d);
fprintf('Pr     %8.2f      %8.2f\n', Pr_u, Pr_d);
fprintf('C/N0   %8.2f      %8.2f\n', CN0_u, CN0_d);
fprintf('LM     %8.2f      %8.2f\n', LM_u, LM_d);

distance = 1000:1000:50000;
Pr = EIRP_d + Gr_d - ...
    (92.45 + 20*log10(f_d) + 20*log10(distance) ...
    + Lrx_d + Latm + Lmisc);

subplot(2,2,1);
plot(distance,Pr,'b','LineWidth',2);
grid on;
title('Received Power vs Distance');
xlabel('Distance (km)');
ylabel('Received Power (dBW)');

E = 20:70;
LM1 = E + Gr_d - LFS_d - Lrx_d - Latm - Lmisc - Pr_req;

subplot(2,2,2);
plot(E,LM1,'r','LineWidth',2);
grid on;
title('Link Margin vs EIRP');
xlabel('EIRP (dBW)');
ylabel('Link Margin (dB)');

L = 180:220;
LM2 = EIRP_d + Gr_d - L - Lrx_d - Latm - Lmisc - Pr_req;

subplot(2,2,3);
plot(L,LM2,'m','LineWidth',2);
grid on;
title('Link Margin vs Free Space Loss');
xlabel('Free Space Loss (dB)');
ylabel('Link Margin (dB)');

G = 10:60;
LM3 = EIRP_d + G - LFS_d - Lrx_d - Latm - Lmisc - Pr_req;

subplot(2,2,4);
plot(G,LM3,'g','LineWidth',2);
%G=10log10(n(piD/lam)^2)
grid on;
title('Link Margin vs Receiver Gain');
xlabel('Receiver Gain (dBi)');
ylabel('Link Margin (dB)');
