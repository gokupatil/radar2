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
