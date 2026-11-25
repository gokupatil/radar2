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
