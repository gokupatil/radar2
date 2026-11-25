function results = exp02_radar_detection_roc(makePlots)
if nargin < 1, makePlots = true; end
savedRng = rng;
rngCleanup = onCleanup(@() rng(savedRng));
rng(2,'twister');
N = 200000;
A = 1;
SNR_dB = [-15 -10 -5 0 5 10 15];
threshold = 0.5;
noiseStd = A./sqrt(10.^(SNR_dB/10));
Pd = zeros(size(SNR_dB));
Pfa = Pd;
for j = 1:numel(SNR_dB)
    r0 = noiseStd(j)*randn(1,N);
    r1 = A + noiseStd(j)*randn(1,N);
    Pfa(j) = mean(r0 > threshold);
    Pd(j) = mean(r1 > threshold);
end
PfaTheory = 0.5*erfc(threshold./(sqrt(2)*noiseStd));
PdTheory = 0.5*erfc((threshold-A)./(sqrt(2)*noiseStd));
assert(all(abs(Pd-PdTheory)<0.015) && all(abs(Pfa-PfaTheory)<0.015));
PfaROC = logspace(-6,0,500);
PdROC = zeros(numel(SNR_dB),numel(PfaROC));
for j = 1:numel(SNR_dB)
    thresholds = sqrt(2)*noiseStd(j)*erfcinv(2*PfaROC);
    PdROC(j,:) = 0.5*erfc((thresholds-A)/(sqrt(2)*noiseStd(j)));
end
summary = table(SNR_dB(:),Pd(:),PdTheory(:),Pfa(:),PfaTheory(:), ...
    'VariableNames',{'SNR_dB','Pd_MC','Pd_Theory','Pfa_MC','Pfa_Theory'});
disp('Experiment 2: real Gaussian amplitude detector, SNR = A^2/noise variance');
disp(summary);
if makePlots
    figure('Name','Experiment 2: Detection Probabilities','Color','w');
    subplot(2,1,1); plot(SNR_dB,Pd,'o',SNR_dB,PdTheory,'-','LineWidth',1.5); grid on;
    xlabel('SNR (dB)'); ylabel('P_d'); legend('Monte Carlo','Theory','Location','best');
    subplot(2,1,2); plot(SNR_dB,Pfa,'o',SNR_dB,PfaTheory,'-','LineWidth',1.5); grid on;
    xlabel('SNR (dB)'); ylabel('P_{fa}'); legend('Monte Carlo','Theory','Location','best');
    figure('Name','Experiment 2: ROC','Color','w');
    semilogx(PfaROC,PdROC.','LineWidth',1.5); grid on; xlim([1e-6 1]); ylim([0 1]);
    xlabel('False-alarm probability'); ylabel('Detection probability'); title('Real Gaussian amplitude-detector ROC');
    legend(arrayfun(@(s) sprintf('%g dB',s),SNR_dB,'UniformOutput',false),'Location','southeast');
    tx = [zeros(1,20) ones(1,20) zeros(1,20)];
    low = tx + 10^(5/20)*randn(size(tx));
    high = tx + 10^(-10/20)*randn(size(tx));
    figure('Name','Experiment 2: Pulse in Noise','Color','w');
    subplot(3,1,1); plot(tx); grid on; title('Transmitted pulse'); ylabel('Amplitude');
    subplot(3,1,2); plot(low); hold on; plot([1 numel(tx)],[threshold threshold],'k--'); grid on; title('Pulse at -5 dB SNR'); ylabel('Amplitude');
    subplot(3,1,3); plot(high); hold on; plot([1 numel(tx)],[threshold threshold],'k--'); grid on; title('Pulse at 10 dB SNR'); xlabel('Sample'); ylabel('Amplitude');
end
results = struct('summary',summary,'PfaROC',PfaROC,'PdROC',PdROC);
end
