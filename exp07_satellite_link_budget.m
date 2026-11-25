function results = exp07_satellite_link_budget(makePlots)
if nargin < 1, makePlots = true; end
frequency_GHz = [14 12];
distance_km = [36000 36000];
Pt_W = [20 10];
Gt_dBi = [45 35];
Gr_dBi = [30 40];
loss_dB = [2 2];
noise_dBW = [-125 -125];
requiredCN_dB = 10;
Pt_dBW = 10*log10(Pt_W);
EIRP = Pt_dBW+Gt_dBi;
FSPL = 92.45+20*log10(frequency_GHz)+20*log10(distance_km);
Pr = EIRP-FSPL+Gr_dBi-loss_dB;
CN = Pr-noise_dBW;
margin = CN-requiredCN_dB;
combinedCN = -10*log10(sum(10.^(-CN/10)));
Parameter = {'Frequency (GHz)';'Distance (km)';'Transmit power (dBW)';'Tx gain (dBi)';'EIRP (dBW)';'Free-space loss (dB)';'Rx gain (dBi)';'Other losses (dB)';'Received power (dBW)';'Noise power (dBW)';'C/N (dB)';'Link margin (dB)'};
values = [frequency_GHz;distance_km;Pt_dBW;Gt_dBi;EIRP;FSPL;Gr_dBi;loss_dB;Pr;noise_dBW;CN;margin];
summary = table(Parameter,values(:,1),values(:,2),'VariableNames',{'Parameter','Uplink','Downlink'});
disp('Experiment 7: satellite link budget'); disp(summary);
fprintf('Combined C/N for independent link noise = %.4f dB\n',combinedCN);
fprintf('Uplink meets requirement: %d; downlink meets requirement: %d\n',margin>=0);
assert(all(isfinite(values(:))) && combinedCN<min(CN));
d = linspace(30000,45000,100);
p = EIRP(2)-(92.45+20*log10(frequency_GHz(2))+20*log10(d))+Gr_dBi(2)-loss_dB(2);
if makePlots
    figure('Name','Experiment 7: Link Budget','Color','w');
    subplot(2,1,1); plot(d,p,'LineWidth',1.5); grid on;
    xlabel('Distance (km)'); ylabel('Received power (dBW)'); title('Downlink received power');
    subplot(2,1,2); bar(margin); set(gca,'XTick',1:2,'XTickLabel',{'Uplink','Downlink'}); grid on;
    ylabel('Margin (dB)'); title('Negative margin means the specified link does not close');
end
results = struct('summary',summary,'CN_dB',CN,'margin_dB',margin,'combinedCN_dB',combinedCN);
end
