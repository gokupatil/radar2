clc; clear; close all;

range = 0:0.1:100;
pulse = @(amp, pos, width) amp * exp(-((range - pos) / width).^2);

clutter    = pulse(8, 20, 2.5) + pulse(6, 75, 3);
statTarget = pulse(7, 50, 1.2);
movTarget1 = pulse(7, 35, 1.2);
movTarget2 = pulse(7, 40, 1.2);

sweep1 = clutter + statTarget + movTarget1;
sweep2 = clutter + statTarget + movTarget2;
mtiOut = abs(sweep1 - sweep2);

figure('Name', 'MTI Radar: Clutter vs Stationary vs Moving', 'Color', 'w');

subplot(3, 1, 1);
plot(range, sweep1, 'b', 'LineWidth', 1.5);
title('1. Sweep 1: Clutter (20, 75 km), Stationary Target (50 km), Moving Target (35 km)');
xlabel('Distance / Range (km)'); ylabel('Amplitude (V)');
ylim([0, 10]); grid on;

subplot(3, 1, 2);
plot(range, sweep2, 'r', 'LineWidth', 1.5);
title('2. Sweep 2 (Delayed): Clutter & Stationary Target Unchanged, Moving Target Shifted to 40 km');
xlabel('Distance / Range (km)'); ylabel('Amplitude (V)');
ylim([0, 10]); grid on;

subplot(3, 1, 3);
plot(range, mtiOut, 'y', 'LineWidth', 1.8);
title('3. MTI Output: Clutter & Stationary Targets Cancelled (0 V), Only Moving Target Detected');
xlabel('Distance / Range (km)'); ylabel('Amplitude (V)');
ylim([0, 10]); grid on;
