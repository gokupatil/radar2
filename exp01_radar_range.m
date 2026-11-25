function results = exp01_radar_range(makePlots)
if nargin < 1, makePlots = true; end
c = 3e8;
Pt = 10e3;
f = 10e9;
Gt = 10^(30/10);
Gr = Gt;
sigma = 1;
L = 10^(6/10);
Pmin = 1e-13;
lambda = c/f;
K = Pt*Gt*Gr*lambda^2*sigma/((4*pi)^3*L);
Rmax = (K/Pmin)^0.25;
R = linspace(1,Rmax,1000);
Pr = K./R.^4;
PtSweep = linspace(1e3,50e3,100);
Rpower = (PtSweep*Gt*Gr*lambda^2*sigma/((4*pi)^3*L*Pmin)).^0.25;
Gain_dBi = linspace(10,50,100);
Rgain = (Pt*10.^(2*Gain_dBi/10)*lambda^2*sigma/((4*pi)^3*L*Pmin)).^0.25;
fSweep = linspace(1,20,100)*1e9;
Rfrequency = (Pt*Gt*Gr*(c./fSweep).^2*sigma/((4*pi)^3*L*Pmin)).^0.25;
sigmaSweep = linspace(0.1,10,100);
Rrcs = (Pt*Gt*Gr*lambda^2*sigmaSweep/((4*pi)^3*L*Pmin)).^0.25;
assert(abs(Pr(end)/Pmin-1) < 1e-12);
assert(all(diff(Pr)<0) && all(diff(Rpower)>0) && all(diff(Rfrequency)<0));
fprintf('Experiment 1: maximum range = %.4f km\n',Rmax/1000);
if makePlots
    figure('Name','Experiment 1: Radar Range','Color','w');
    subplot(3,2,1); plot(R/1000,10*log10(Pr/1e-3),'LineWidth',1.5); grid on;
    xlabel('Range (km)'); ylabel('Received power (dBm)'); title('Received power vs range');
    subplot(3,2,2); plot(PtSweep/1000,Rpower/1000,'LineWidth',1.5); grid on;
    xlabel('Transmit power (kW)'); ylabel('Maximum range (km)'); title('Transmit power sweep');
    subplot(3,2,3); plot(Gain_dBi,Rgain/1000,'LineWidth',1.5); grid on;
    xlabel('Tx and Rx gain (dBi)'); ylabel('Maximum range (km)'); title('Equal antenna gains');
    subplot(3,2,4); plot(fSweep/1e9,Rfrequency/1000,'LineWidth',1.5); grid on;
    xlabel('Frequency (GHz)'); ylabel('Maximum range (km)'); title('Fixed antenna gains');
    subplot(3,2,5); plot(sigmaSweep,Rrcs/1000,'LineWidth',1.5); grid on;
    xlabel('Radar cross section (m^2)'); ylabel('Maximum range (km)'); title('RCS sweep');
end
results = struct('maximumRange_m',Rmax,'range_m',R,'receivedPower_W',Pr);
end
