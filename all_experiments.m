function results = all_experiments(makePlots,animate)
if nargin < 1, makePlots = true; end
if nargin < 2, animate = true; end
results = cell(1,9);
results{1} = exp01_radar_range(makePlots);
results{2} = exp02_radar_detection_roc(makePlots);
results{3} = exp03_doppler_shift(makePlots);
results{4} = exp04_mti_canceller(makePlots);
results{5} = exp05_pulse_doppler(makePlots);
results{6} = exp06_radar_display(makePlots,animate);
results{7} = exp07_satellite_link_budget(makePlots);
results{8} = exp08_satellite_orbits(makePlots,animate);
results{9} = exp09_system_noise(makePlots);
disp('All nine experiments completed and internal checks passed.');
end

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

function results = exp03_doppler_shift(makePlots)
if nargin < 1, makePlots = true; end
savedRng = rng;
rngCleanup = onCleanup(@() rng(savedRng));
rng(3,'twister');
c = 3e8;
fc = 10e9;
lambda = c/fc;
velocity = 0:10:100;
fdVelocity = 2*velocity/lambda;
fcSweep = [1 3 5 10 24]*1e9;
fdFrequency = 2*50*fcSweep/c;
fs = 1e5;
N = 1000;
t = (0:N-1)/fs;
fd = 2*60/lambda;
baseband = exp(1i*2*pi*fd*t) + 0.5/sqrt(2)*(randn(1,N)+1i*randn(1,N));
frequency = (-N/2:N/2-1)*fs/N;
PSD = abs(fftshift(fft(baseband))).^2/(fs*N);
[~,peak] = max(PSD);
estimate = frequency(peak);
assert(abs(estimate-fd) <= fs/N);
disp('Experiment 3: positive velocity denotes positive Doppler');
disp(table(velocity(:),fdVelocity(:),'VariableNames',{'Velocity_m_s','Doppler_Hz'}));
disp(table(fcSweep(:)/1e9,fdFrequency(:),'VariableNames',{'Carrier_GHz','Doppler_Hz'}));
fprintf('Expected Doppler = %.2f Hz; FFT estimate = %.2f Hz\n',fd,estimate);
if makePlots
    figure('Name','Experiment 3: Doppler Shift','Color','w');
    subplot(3,1,1); plot(velocity,fdVelocity,'-o'); grid on; xlabel('Radial velocity (m/s)'); ylabel('Doppler (Hz)');
    subplot(3,1,2); plot(fcSweep/1e9,fdFrequency,'-o'); grid on; xlabel('Carrier frequency (GHz)'); ylabel('Doppler (Hz)');
    subplot(3,1,3); plot(frequency,10*log10(max(PSD,realmin))); grid on; xlim([-10000 10000]);
    xlabel('Baseband frequency (Hz)'); ylabel('PSD (dB/Hz)'); title('Complex baseband Doppler spectrum');
end
results = struct('expectedDoppler_Hz',fd,'estimatedDoppler_Hz',estimate,'frequency_Hz',frequency,'PSD',PSD);
end

function results = exp04_mti_canceller(makePlots)
if nargin < 1, makePlots = true; end
savedRng = rng;
rngCleanup = onCleanup(@() rng(savedRng));
rng(4,'twister');
PRF = 1000;
N = 200;
t = (0:N-1)/PRF;
fd = 100;
clutter = 2*ones(1,N);
target = 4*sin(2*pi*fd*t);
rx = clutter + target + 0.05*randn(1,N);
y = diff(rx);
clutterResidual = diff(clutter);
assert(all(clutterResidual==0));
Nfft = 1024;
f = (-Nfft/2:Nfft/2-1)*PRF/Nfft;
X = abs(fftshift(fft(rx(2:end),Nfft)))/(N-1);
Y = abs(fftshift(fft(y,Nfft)))/(N-1);
fResponse = linspace(-PRF,PRF,2001);
H = abs(1-exp(-1i*2*pi*fResponse/PRF));
gain = 2*abs(sin(pi*fd/PRF));
fprintf('Experiment 4: clutter residual = %g; target amplitude gain = %.6f\n',max(abs(clutterResidual)),gain);
if makePlots
    figure('Name','Experiment 4: MTI','Color','w');
    subplot(3,1,1); plot(t,rx); grid on; xlabel('Slow time (s)'); ylabel('Amplitude'); title('Received clutter, target and noise');
    subplot(3,1,2); plot(t,clutter,t,target); grid on; legend('Clutter','Moving target'); xlabel('Slow time (s)'); ylabel('Amplitude');
    subplot(3,1,3); plot(t(2:end),y); grid on; xlabel('Slow time (s)'); ylabel('Amplitude'); title('Two-pulse canceller output');
    figure('Name','Experiment 4: MTI Spectra','Color','w');
    subplot(2,1,1); plot(f,X,f,Y); grid on; xlabel('Doppler frequency (Hz)'); ylabel('Magnitude'); legend('Before MTI','After MTI');
    subplot(2,1,2); plot(fResponse,H); grid on; xlabel('Doppler frequency (Hz)'); ylabel('|H|'); title('Nulls at integer multiples of PRF');
end
results = struct('output',y,'clutterResidual',clutterResidual,'targetGain',gain,'PRF_Hz',PRF);
end

function results = exp05_pulse_doppler(makePlots)
if nargin < 1, makePlots = true; end
savedRng = rng;
rngCleanup = onCleanup(@() rng(savedRng));
rng(5,'twister');
c = 3e8;
fc = 10e9;
lambda = c/fc;
PRF = 2000;
Fs = 20e6;
pulseWidth = 1e-6;
bandwidth = 5e6;
Np = 128;
Ns = round(Fs/PRF);
L = round(pulseWidth*Fs);
fastTime = (0:L-1)'/Fs;
tx = exp(1i*pi*(bandwidth/pulseWidth)*(fastTime-pulseWidth/2).^2);
targetRanges = [3000 6000];
targetVelocities = [12 -8];
assert(all(abs(targetVelocities)<lambda*PRF/4));
rx = complex(zeros(Ns,Np));
slowTime = (0:Np-1)/PRF;
for j = 1:numel(targetRanges)
    delay = round(2*targetRanges(j)*Fs/c);
    assert(delay+L <= Ns);
    fd = 2*targetVelocities(j)/lambda;
    echo = tx.*exp(1i*2*pi*fd*fastTime);
    rx(delay+(1:L),:) = rx(delay+(1:L),:) + echo*exp(1i*2*pi*fd*slowTime);
end
rx = rx + 0.1*(randn(Ns,Np)+1i*randn(Ns,Np));
Nfft = 2^nextpow2(Ns+L-1);
h = conj(flipud(tx));
matchedFull = ifft(bsxfun(@times,fft(rx,Nfft,1),fft(h,Nfft)),[],1);
matched = matchedFull(L:L+Ns-1,:);
window = 0.5-0.5*cos(2*pi*(0:Np-1)/(Np-1));
RD = fftshift(fft(bsxfun(@times,matched,window),[],2),2);
rangeAxis = (0:Ns-1)'*c/(2*Fs);
fdAxis = (-Np/2:Np/2-1)*PRF/Np;
velocityAxis = fdAxis*lambda/2;
map_dB = 20*log10(max(abs(RD)/max(abs(RD(:))),1e-6));
estimatedRanges = zeros(size(targetRanges));
estimatedVelocities = zeros(size(targetRanges));
for j = 1:numel(targetRanges)
    rows = find(abs(rangeAxis-targetRanges(j))<=2*c/(2*bandwidth));
    block = abs(RD(rows,:));
    [~,peak] = max(block(:));
    [ir,iv] = ind2sub(size(block),peak);
    estimatedRanges(j) = rangeAxis(rows(ir));
    estimatedVelocities(j) = velocityAxis(iv);
end
assert(all(abs(estimatedRanges-targetRanges)<=c/(2*bandwidth)));
assert(all(abs(estimatedVelocities-targetVelocities)<=lambda*PRF/(2*Np)));
summary = table(targetRanges(:),estimatedRanges(:),targetVelocities(:),estimatedVelocities(:), ...
    'VariableNames',{'TrueRange_m','EstimatedRange_m','TrueVelocity_m_s','EstimatedVelocity_m_s'});
disp('Experiment 5: range-Doppler target peaks'); disp(summary);
fprintf('Range resolution = %.2f m; range sampling = %.2f m; velocity bin = %.6f m/s\n',c/(2*bandwidth),c/(2*Fs),lambda*PRF/(2*Np));
if makePlots
    figure('Name','Experiment 5: Range-Doppler Map','Color','w');
    imagesc(velocityAxis,rangeAxis,map_dB); axis xy; ylim([0 8000]); caxis([-60 0]); colorbar;
    xlabel('Radial velocity (m/s)'); ylabel('Range (m)'); title('Range-Doppler magnitude (dB relative to peak)');
end
results = struct('summary',summary,'range_m',rangeAxis,'velocity_m_s',velocityAxis,'map_dB',map_dB);
end

function results = exp06_radar_display(makePlots,animate)
if nargin < 1, makePlots = true; end
if nargin < 2, animate = true; end
savedRng = rng;
rngCleanup = onCleanup(@() rng(savedRng));
rng(6,'twister');
ranges = [30 75 50];
azimuth = [45 120 300];
amplitude = [0.8 0.5 0.9];
vr = [-0.05 -0.08 0.04];
vaz = [0.1 -0.15 0.2];
maxRange = 100;
step = 2;
numFrames = 5*360/step;
x = ranges.*cosd(azimuth);
y = ranges.*sind(azimuth);
initial = table((1:3)',ranges',azimuth',x',y', ...
    'VariableNames',{'Target','Range_km','Azimuth_deg','X_km','Y_km'});
disp('Experiment 6: azimuth is counterclockwise from +X; motion is per frame'); disp(initial);
rangeGrid = 0:0.1:maxRange;
detectionCount = zeros(1,3);
timer = zeros(1,3);
if makePlots
    fig = figure('Name','Experiment 6: Radar Displays','Color','w');
    ax1 = subplot(2,2,1); aLine = plot(rangeGrid,zeros(size(rangeGrid))); grid on;
    xlim([0 maxRange]); ylim([-0.2 1.1]); xlabel('Range (km)'); ylabel('Amplitude'); title('A-scan at current beam angle');
    ax2 = subplot(2,2,3); bPoints = scatter(azimuth,ranges,60,amplitude,'filled'); grid on;
    xlim([0 360]); ylim([0 maxRange]); xlabel('Azimuth (degrees)'); ylabel('Range (km)'); title('B-scan: simulated target positions');
    ax3 = subplot(2,2,[2 4]); hold on; axis equal; grid on;
    theta = linspace(0,2*pi,361);
    for radius = 20:20:maxRange
        plot(radius*cos(theta),radius*sin(theta),'Color',[0.75 0.75 0.75]);
    end
    xlim([-maxRange maxRange]); ylim([-maxRange maxRange]); xlabel('X (km)'); ylabel('Y (km)');
    title('PPI: red truth, green recent detections');
    sweepLine = plot([0 maxRange],[0 0],'g-','LineWidth',1.5);
    truthPoints = scatter(x,y,35,'r','filled');
    hitPoints = scatter(nan(1,3),nan(1,3),100,'g');
end
for frame = 1:numFrames
    sweep = mod((frame-1)*step,360);
    ranges = ranges + vr;
    above = ranges > maxRange;
    below = ranges < 0;
    ranges(above) = 2*maxRange-ranges(above);
    ranges(below) = -ranges(below);
    vr(above | below) = -vr(above | below);
    azimuth = mod(azimuth+vaz,360);
    separation = abs(mod(azimuth-sweep+180,360)-180);
    hits = separation <= step;
    detectionCount = detectionCount + hits;
    timer = max(timer-1,0);
    timer(hits) = 25;
    x = ranges.*cosd(azimuth);
    y = ranges.*sind(azimuth);
    ascan = 0.025*randn(size(rangeGrid));
    for j = find(hits)
        ascan = ascan + amplitude(j)*exp(-0.5*((rangeGrid-ranges(j))/0.2).^2);
    end
    if makePlots && (animate || frame==numFrames)
        if ~ishghandle(fig), break; end
        set(aLine,'YData',ascan);
        set(bPoints,'XData',azimuth,'YData',ranges);
        set(sweepLine,'XData',[0 maxRange*cosd(sweep)],'YData',[0 maxRange*sind(sweep)]);
        set(truthPoints,'XData',x,'YData',y);
        hx = x; hy = y;
        hx(timer==0) = NaN; hy(timer==0) = NaN;
        set(hitPoints,'XData',hx,'YData',hy);
        title(ax1,sprintf('A-scan at %.0f degrees',sweep));
        drawnow;
        if animate, pause(0.005); end
    end
end
assert(all(ranges>=0 & ranges<=maxRange));
assert(max(abs(hypot(x,y)-ranges))<1e-10);
results = struct('initial',initial,'finalRange_km',ranges,'finalAzimuth_deg',azimuth,'detectionCount',detectionCount,'framesCompleted',frame);
end

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

function results = exp08_satellite_orbits(makePlots,animate)
if nargin < 1, makePlots = true; end
if nargin < 2, animate = true; end
aEarth = 6378.137;
bEarth = aEarth*(1-1/298.257223563);
mu = 398600.4418;
omegaEarth = 7.2921150e-5;
names = {'LEO','MEO','GEO'};
altitude = [500 20200 35786];
a = aEarth+altitude;
ecc = [0.001 0.010 0.0001];
inc = [51.6 55 0]*pi/180;
raan = [0 0 0];
argPeri = [0 0 0];
M0 = [0 45 0]*pi/180;
n = sqrt(mu./a.^3);
period = 2*pi./n;
summary = table(names',altitude',a',period'/3600,ecc', ...
    'VariableNames',{'Orbit','NominalAltitude_km','SemimajorAxis_km','Period_hours','Eccentricity'});
disp('Experiment 8: two-body Keplerian orbits; Earth rotation angle is zero at t = 0'); disp(summary);
time = linspace(0,24*3600,3000);
tracks = cell(1,3);
colors = lines(3);
if makePlots
    orbitFig = figure('Name','Experiment 8: Inertial Orbits','Color','w');
    orbitAx = axes; hold(orbitAx,'on');
    [sx,sy,sz] = sphere(50);
    surf(orbitAx,aEarth*sx,aEarth*sy,bEarth*sz,'FaceColor',[0.15 0.4 0.7],'EdgeColor','none','HandleVisibility','off');
    satellites = gobjects(1,3);
end
for j = 1:3
    r = orbit_position(a(j),ecc(j),inc(j),raan(j),argPeri(j),M0(j)+n(j)*time);
    angle = omegaEarth*time;
    x = cos(angle).*r(1,:) + sin(angle).*r(2,:);
    y = -sin(angle).*r(1,:) + cos(angle).*r(2,:);
    z = r(3,:);
    lon = atan2d(y,x);
    lat = orbit_latitude(x,y,z,aEarth,bEarth);
    radius = sqrt(sum(r.^2,1));
    assert(all(radius>=a(j)*(1-ecc(j))-1e-7) && all(radius<=a(j)*(1+ecc(j))+1e-7));
    assert(all(abs(lat)<=90) && all(abs(lon)<=180));
    tracks{j} = struct('longitude_deg',lon,'latitude_deg',lat,'radius_km',radius);
    if makePlots
        path = orbit_position(a(j),ecc(j),inc(j),raan(j),argPeri(j),linspace(0,2*pi,700));
        plot3(orbitAx,path(1,:),path(2,:),path(3,:),'Color',colors(j,:),'DisplayName',names{j});
        satellites(j) = plot3(orbitAx,r(1,1),r(2,1),r(3,1),'o','MarkerFaceColor',colors(j,:),'Color',colors(j,:),'HandleVisibility','off');
    end
end
if makePlots
    axis(orbitAx,'equal'); grid(orbitAx,'on'); view(orbitAx,35,24);
    xlabel(orbitAx,'Inertial X (km)'); ylabel(orbitAx,'Inertial Y (km)'); zlabel(orbitAx,'Inertial Z (km)');
    title(orbitAx,'Keplerian orbits'); legend(orbitAx,'Location','best');
    figure('Name','Experiment 8: Ground Tracks','Color','w');
    for j = 1:3
        lon = tracks{j}.longitude_deg; lat = tracks{j}.latitude_deg;
        seam = find(abs(diff(lon))>180)+1;
        lon(seam) = NaN; lat(seam) = NaN;
        subplot(3,1,j); plot(lon,lat,'Color',colors(j,:)); grid on;
        xlim([-180 180]); ylim([-90 90]); xlabel('Longitude (degrees)'); ylabel('Geodetic latitude (degrees)');
        title([names{j} ' ground track over 24 hours']);
    end
    if animate
        for t = linspace(0,24*3600,360)
            if ~ishghandle(orbitFig), break; end
            for j = 1:3
                r = orbit_position(a(j),ecc(j),inc(j),raan(j),argPeri(j),M0(j)+n(j)*t);
                set(satellites(j),'XData',r(1),'YData',r(2),'ZData',r(3));
            end
            drawnow; pause(0.005);
        end
    end
end
results = struct('summary',summary,'time_s',time,'groundTracks',{tracks});
end

function r = orbit_position(a,e,i,raan,argPeri,M)
M = mod(M+pi,2*pi)-pi;
E = M;
for iteration = 1:30
    delta = (E-e*sin(E)-M)./(1-e*cos(E));
    E = E-delta;
    if max(abs(delta(:)))<1e-13, break; end
end
assert(max(abs(E-e*sin(E)-M))<1e-11);
rp = [a*(cos(E)-e);a*sqrt(1-e^2)*sin(E);zeros(size(E))];
Rz1 = [cos(raan) -sin(raan) 0;sin(raan) cos(raan) 0;0 0 1];
Rx = [1 0 0;0 cos(i) -sin(i);0 sin(i) cos(i)];
Rz2 = [cos(argPeri) -sin(argPeri) 0;sin(argPeri) cos(argPeri) 0;0 0 1];
r = Rz1*Rx*Rz2*rp;
end

function latitude = orbit_latitude(x,y,z,a,b)
e2 = 1-b^2/a^2;
p = hypot(x,y);
lat = atan2(z,p*(1-e2));
for iteration = 1:12
    N = a./sqrt(1-e2*sin(lat).^2);
    next = atan2(z+e2*N.*sin(lat),p);
    change = max(abs(next(:)-lat(:)));
    lat = next;
    if change<1e-13, break; end
end
latitude = lat*180/pi;
end

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
