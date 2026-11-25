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
