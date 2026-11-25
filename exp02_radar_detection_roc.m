clc; clear; close all;

snr = [5 7.5 10 12.5];
T = linspace(0,6,50);
x = zeros(1,500);
x(100:150) = 1;

N = 1e5;
Pfa = zeros(4,length(T));
Pd = zeros(4,length(T));

figure('Name','Noisy vs Clean Signals','Color','w');

for k = 1:4
    subplot(4,1,k);
    plot(awgn(x,snr(k),'measured')); hold on;
    plot(x,'r','LineWidth',1.5);
    grid on;
    title("SNR = " + snr(k) + " dB");
    xlabel('Sample');
    ylabel('Amplitude');
    if k == 1
        legend('Noisy','Clean');
    end
end

for k = 1:4
    s = 10^(-snr(k)/20);
    Pfa(k,:) = mean(s*randn(N,1) > T,1);
    Pd(k,:) = mean(1 + s*randn(N,1) > T,1);
end

figure('Name','Pd and Pfa vs Threshold','Color','w');
plot(T,Pd',T,Pfa','--');
grid on;
xlabel('Threshold (T)');
ylabel('Probability');
title('Pd and Pfa vs Threshold');
legend(["Pd " + snr + " dB", "Pfa " + snr + " dB"]);

figure('Name','ROC Curves','Color','w');
plot(Pfa',Pd','-o');
grid on;
xlabel('Probability of False Alarm (Pfa)');
ylabel('Probability of Detection (Pd)');
title('ROC Curves');
legend("SNR = " + snr + " dB",'Location','southeast');

fixedThreshold = 0.5;
Pd_SNR = zeros(1,4);
Pfa_SNR = zeros(1,4);

for k = 1:4
    s = 10^(-snr(k)/20);
    noise = s*randn(N,1);

    Pfa_SNR(k) = mean(noise > fixedThreshold);
    Pd_SNR(k) = mean(1 + noise > fixedThreshold);
end

figure('Name','SNR vs Pd and Pfa','Color','w');

plot(snr,Pd_SNR,'b-o','LineWidth',1.5,...
    'MarkerFaceColor','b');
hold on;
plot(snr,Pfa_SNR,'r-s','LineWidth',1.5,...
    'MarkerFaceColor','r');

grid on;
xlabel('SNR (dB)');
ylabel('Probability');
title('SNR vs Pd and Pfa');
legend('Probability of Detection (Pd)',...
    'Probability of False Alarm (Pfa)',...
    'Location','best');

ylim([0 1]);
