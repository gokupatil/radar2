clc;clear;close all;

rng(5);

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

t = (0:L-1)'/Fs;
tx = exp(1i*pi*(bandwidth/pulseWidth)*(t-pulseWidth/2).^2);

targetRanges = [3000 6000];
targetVelocities = [12 -8];

rx = complex(zeros(Ns,Np));
slowTime = (0:Np-1)/PRF;

for j = 1:2
    delay = round(2*targetRanges(j)*Fs/c);
    fd = 2*targetVelocities(j)/lambda;

    echo = tx.*exp(1i*2*pi*fd*t);

    rx(delay+(1:L),:) = rx(delay+(1:L),:) + ...
        echo*exp(1i*2*pi*fd*slowTime);
end

rx = rx + 0.1*(randn(Ns,Np)+1i*randn(Ns,Np));

h = conj(flipud(tx));
Nfft = 2^nextpow2(Ns+L-1);

matchedFull = ifft(fft(rx,Nfft,1).*fft(h,Nfft),[],1);
matched = matchedFull(L:L+Ns-1,:);

window = 0.5-0.5*cos(2*pi*(0:Np-1)/(Np-1));
RD = fftshift(fft(matched.*window,[],2),2);

rangeAxis = (0:Ns-1)'*c/(2*Fs);
fdAxis = (-Np/2:Np/2-1)*PRF/Np;
velocityAxis = fdAxis*lambda/2;

map_dB = 20*log10(max(abs(RD)/max(abs(RD(:))),1e-6));

estimatedRanges = zeros(1,2);
estimatedVelocities = zeros(1,2);

for j = 1:2
    rows = find(abs(rangeAxis-targetRanges(j)) <= 2*c/(2*bandwidth));
    block = abs(RD(rows,:));

    [~,peak] = max(block(:));
    [ir,iv] = ind2sub(size(block),peak);

    estimatedRanges(j) = rangeAxis(rows(ir));
    estimatedVelocities(j) = velocityAxis(iv);
end

summary = table(targetRanges(:),estimatedRanges(:), ...
    targetVelocities(:),estimatedVelocities(:), ...
    'VariableNames',{'TrueRange_m','EstimatedRange_m', ...
    'TrueVelocity_m_s','EstimatedVelocity_m_s'});

disp('Experiment 5: Range-Doppler Target Peaks');
disp(summary);

fprintf('Range resolution = %.2f m\n',c/(2*bandwidth));
fprintf('Range sampling = %.2f m\n',c/(2*Fs));
fprintf('Velocity bin = %.6f m/s\n',lambda*PRF/(2*Np));

figure('Name','Experiment 5: Range-Doppler Map','Color','w');

imagesc(velocityAxis,rangeAxis,map_dB);
axis xy;
ylim([0 8000]);
caxis([-60 0]);
colorbar;

xlabel('Radial velocity (m/s)');
ylabel('Range (m)');
title('Range-Doppler Magnitude (dB relative to peak)');
