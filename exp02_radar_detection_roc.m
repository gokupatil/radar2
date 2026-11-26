clc; clear; close all;

snr = [5 7.5 10 12.5];
s = 10.^(-snr/20);
T = linspace(0,6,50);
N = 1e5;

x = zeros(1,500);
x(100:150) = 1;

figure('Name','Noisy vs Clean Signals','Color','w');
subplot(4,1,1); plot(awgn(x,snr(1),'measured')); hold on; plot(x,'r','LineWidth',1.5); grid on; title('SNR = 5 dB'); xlabel('Sample'); ylabel('Amplitude'); legend('Noisy','Clean');
subplot(4,1,2); plot(awgn(x,snr(2),'measured')); hold on; plot(x,'r','LineWidth',1.5); grid on; title('SNR = 7.5 dB'); xlabel('Sample'); ylabel('Amplitude');
subplot(4,1,3); plot(awgn(x,snr(3),'measured')); hold on; plot(x,'r','LineWidth',1.5); grid on; title('SNR = 10 dB'); xlabel('Sample'); ylabel('Amplitude');
subplot(4,1,4); plot(awgn(x,snr(4),'measured')); hold on; plot(x,'r','LineWidth',1.5); grid on; title('SNR = 12.5 dB'); xlabel('Sample'); ylabel('Amplitude');

noise = randn(N,1) .* s;
noise3 = permute(noise,[1 3 2]);
Pfa = squeeze(mean(noise3 > T,1))';
Pd = squeeze(mean(1 + noise3 > T,1))';

figure('Name','Pd and Pfa vs Threshold','Color','w');
plot(T,Pd',T,Pfa','--');
grid on;
xlabel('Threshold (T)'); ylabel('Probability');
title('Pd and Pfa vs Threshold');
legend(["Pd " + snr + " dB", "Pfa " + snr + " dB"]);

figure('Name','ROC Curves','Color','w');
plot(Pfa',Pd','-o');
grid on;
xlabel('Probability of False Alarm (Pfa)'); ylabel('Probability of Detection (Pd)');
title('ROC Curves');
legend("SNR = " + snr + " dB",'Location','southeast');

Pfa_SNR = mean(noise > 0.5);
Pd_SNR = mean(1 + noise > 0.5);

figure('Name','SNR vs Pd and Pfa','Color','w');
plot(snr,Pd_SNR,'b-o','LineWidth',1.5,'MarkerFaceColor','b'); hold on;
plot(snr,Pfa_SNR,'r-s','LineWidth',1.5,'MarkerFaceColor','r');
grid on;
xlabel('SNR (dB)'); ylabel('Probability');
title('SNR vs Pd and Pfa');
legend('Probability of Detection (Pd)','Probability of False Alarm (Pfa)','Location','best');
ylim([0 1]);
