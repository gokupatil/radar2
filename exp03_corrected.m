clc; clear; close all;

c  = 3e8;
fc = 10e9;

vel = [10 20 30 40 50];
fdVel = 2 * vel * fc / c;

figure;
plot(vel, fdVel, '-o', 'LineWidth', 2, 'MarkerSize', 7);
grid on;
xlabel('Target Velocity (m/s)');
ylabel('Doppler Frequency (Hz)');
title('Target Velocity vs Doppler Frequency');

opFreq = [2 5 10 15 20] * 1e9;
fdFreq = 2 * 30 * opFreq / c;

figure;
plot(opFreq/1e9, fdFreq, '-o', 'LineWidth', 2, 'MarkerSize', 7);
grid on;
xlabel('Operating Frequency (GHz)');
ylabel('Doppler Frequency (Hz)');
title('Operating Frequency vs Doppler Frequency');

dist = [50 100 200 500];
relPow = 1 ./ dist.^4;
relPow = relPow / max(relPow);

figure;
semilogy(dist, relPow, '-o', 'LineWidth', 2, 'MarkerSize', 7);
grid on;
xlabel('Target Distance (m)');
ylabel('Normalized Received Power');
title('Effect of Target Distance on Received Power');

v = 30; R = 100;
fd = 2 * v * fc / c;

fprintf('Target velocity = %.2f m/s\n', v);
fprintf('Target distance = %.2f m\n', R);
fprintf('Carrier frequency = %.2f GHz\n', fc/1e9);
fprintf('Doppler frequency = %.2f Hz\n', fd);

Fs = 100e3;
T  = 0.1;
t  = 0:1/Fs:T-1/Fs;
N  = length(t);

rxSig    = exp(1j * 2 * pi * fd * t);
noiseVar = mean(abs(rxSig).^2) / 10^(15/10);
rxSig    = rxSig + sqrt(noiseVar/2) * (randn(size(rxSig)) + 1j*randn(size(rxSig)));

f = (-N/2 : N/2-1) * (Fs/N);
P = abs(fftshift(fft(rxSig))).^2 / N;

figure;
subplot(2,1,1);
plot(t, real(rxSig), 'b');
grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Received Doppler Signal');

subplot(2,1,2);
plot(f, 10*log10(P + eps), 'LineWidth', 1.5);
grid on; hold on;
xline(fd, '--r', 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title('Doppler Power Spectrum Using FFT');
xlim([-5000 5000]);
legend('FFT Power Spectrum', sprintf('Expected f_D = %.0f Hz', fd));

[~, idx] = max(P);
detFd = abs(f(idx));

fprintf('\nExpected Doppler frequency = %.2f Hz\n', fd);
fprintf('Detected FFT frequency   = %.2f Hz\n', detFd);
