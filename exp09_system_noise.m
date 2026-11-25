function results = exp09_system_noise(makePlots)
if nargin < 1, makePlots = true; end
k = 1.380649e-23;
C_dBW = -100;
B = 1e6;
T = 150;
G_dBi = 35;
noise_W = k*T*B;
noise_dBW = 10*log10(noise_W);
CN_dB = C_dBW-noise_dBW;
GT_dBK = G_dBi-10*log10(T);
CN0_dBHz = CN_dB+10*log10(B);
temperatures = 50:10:500;
gains = 10:60;
noiseSweep = 10*log10(k*temperatures*B);
cnSweep = C_dBW-noiseSweep;
gtGain = gains-10*log10(T);
gtTemperature = G_dBi-10*log10(temperatures);
assert(all(diff(noiseSweep)>0) && all(diff(cnSweep)<0) && all(diff(gtTemperature)<0));
fprintf('Experiment 9:\nNoise = %.8g W = %.6f dBW\nC/N = %.6f dB\nG/T = %.6f dB/K\nC/N0 = %.6f dB-Hz\n',noise_W,noise_dBW,CN_dB,GT_dBK,CN0_dBHz);
if makePlots
    figure('Name','Experiment 9: System Noise','Color','w');
    subplot(2,2,1); plot(temperatures,noiseSweep,'LineWidth',1.5); grid on; xlabel('System temperature (K)'); ylabel('Noise (dBW)');
    subplot(2,2,2); plot(temperatures,cnSweep,'LineWidth',1.5); grid on; xlabel('System temperature (K)'); ylabel('C/N (dB)');
    subplot(2,2,3); plot(gains,gtGain,'LineWidth',1.5); grid on; xlabel('Antenna gain (dBi)'); ylabel('G/T (dB/K)');
    subplot(2,2,4); plot(temperatures,gtTemperature,'LineWidth',1.5); grid on; xlabel('System temperature (K)'); ylabel('G/T (dB/K)');
end
results = struct('noise_W',noise_W,'noise_dBW',noise_dBW,'CN_dB',CN_dB,'GT_dBK',GT_dBK,'CN0_dBHz',CN0_dBHz);
end
