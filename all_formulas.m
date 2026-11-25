function F = all_formulas()
F.c_m_s = 3e8;
F.k_J_K = 1.380649e-23;
F.mu_km3_s2 = 398600.4418;
F.earthRotation_rad_s = 7.2921150e-5;
c = F.c_m_s;
k = F.k_J_K;
mu = F.mu_km3_s2;
omegaEarth = F.earthRotation_rad_s;
F.wavelength_m = @(frequency_Hz) c./frequency_Hz;
F.dB_to_power_ratio = @(dB) 10.^(dB/10);
F.power_ratio_to_dB = @(ratio) 10*log10(ratio);
F.watts_to_dBW = @(W) 10*log10(W);
F.watts_to_dBm = @(W) 10*log10(W/1e-3);
F.radar_received_W = @(Pt,Gt,Gr,lambda,sigma,R,L) Pt.*Gt.*Gr.*lambda.^2.*sigma./((4*pi)^3.*R.^4.*L);
F.radar_maximum_range_m = @(Pt,Gt,Gr,lambda,sigma,L,Pmin) (Pt.*Gt.*Gr.*lambda.^2.*sigma./((4*pi)^3.*L.*Pmin)).^0.25;
F.real_noise_std = @(A,SNR_dB) abs(A)./sqrt(10.^(SNR_dB/10));
F.real_detector_Pfa = @(threshold,noiseStd) 0.5*erfc(threshold./(sqrt(2)*noiseStd));
F.real_detector_Pd = @(threshold,A,noiseStd) 0.5*erfc((threshold-A)./(sqrt(2)*noiseStd));
F.real_detector_threshold = @(Pfa,noiseStd) sqrt(2)*noiseStd.*erfcinv(2*Pfa);
F.real_detector_ROC = @(Pfa,SNR_dB) 0.5*erfc(erfcinv(2*Pfa)-sqrt(10.^(SNR_dB/10)/2));
F.doppler_Hz = @(velocity_m_s,lambda_m) 2*velocity_m_s./lambda_m;
F.velocity_m_s = @(doppler_Hz,lambda_m) doppler_Hz.*lambda_m/2;
F.PRI_s = @(PRF_Hz) 1./PRF_Hz;
F.two_pulse_canceller = @(x) diff(x,1,2);
F.MTI_response = @(fd,PRF) 1-exp(-1i*2*pi*fd./PRF);
F.MTI_magnitude = @(fd,PRF) 2*abs(sin(pi*fd./PRF));
F.blind_speed_m_s = @(order,lambda,PRF) order.*lambda.*PRF/2;
F.round_trip_delay_s = @(range_m) 2*range_m/c;
F.range_from_delay_m = @(delay_s) c*delay_s/2;
F.range_resolution_m = @(bandwidth_Hz) c./(2*bandwidth_Hz);
F.range_sample_spacing_m = @(sampleRate_Hz) c./(2*sampleRate_Hz);
F.unambiguous_range_m = @(PRF_Hz) c./(2*PRF_Hz);
F.unambiguous_speed_m_s = @(lambda,PRF) lambda.*PRF/4;
F.doppler_bin_Hz = @(PRF,numPulses) PRF./numPulses;
F.velocity_bin_m_s = @(lambda,PRF,numPulses) lambda.*PRF./(2*numPulses);
F.centered_fft_axis_Hz = @(N,Fs) (-floor(N/2):ceil(N/2)-1)*Fs/N;
F.chirp_baseband = @(t,B,T) exp(1i*pi*(B/T)*(t-T/2).^2);
F.matched_filter = @(tx) conj(flipud(tx(:)));
F.linear_matched_output = @(rx,tx) conv(rx(:),conj(flipud(tx(:))));
F.doppler_transform = @(rangePulseMatrix) fftshift(fft(rangePulseMatrix,[],2),2);
F.PSD_two_sided = @(x,Fs) abs(fftshift(fft(x))).^2/(Fs*numel(x));
F.relative_magnitude_dB = @(A) 20*log10(max(abs(A)/max(max(abs(A(:))),realmin),1e-6));
F.polar_x = @(R,azimuth_deg) R.*cosd(azimuth_deg);
F.polar_y = @(R,azimuth_deg) R.*sind(azimuth_deg);
F.circular_separation_deg = @(a,b) abs(mod(a-b+180,360)-180);
F.EIRP_dBW = @(Pt_W,Gt_dBi,txLoss_dB) 10*log10(Pt_W)+Gt_dBi-txLoss_dB;
F.FSPL_dB = @(frequency_GHz,distance_km) 92.45+20*log10(frequency_GHz)+20*log10(distance_km);
F.link_received_dBW = @(EIRP,FSPL,Gr,otherLoss) EIRP-FSPL+Gr-otherLoss;
F.CN_dB = @(carrier_dBW,noise_dBW) carrier_dBW-noise_dBW;
F.margin_dB = @(CN,requiredCN) CN-requiredCN;
F.combined_CN_dB = @(CN_dB) -10*log10(sum(10.^(-CN_dB(:)/10)));
F.mean_motion_rad_s = @(a_km) sqrt(mu./a_km.^3);
F.orbit_period_s = @(a_km) 2*pi*sqrt(a_km.^3/mu);
F.mean_anomaly_rad = @(M0,n,t) M0+n.*t;
F.kepler_residual = @(E,e,M) E-e.*sin(E)-M;
F.orbit_radius_km = @(a,e,E) a.*(1-e.*cos(E));
F.orbit_ECI_km = @orbit_position;
F.ECI_to_ECEF_km = @(r,t) [cos(omegaEarth*t).*r(1,:)+sin(omegaEarth*t).*r(2,:);-sin(omegaEarth*t).*r(1,:)+cos(omegaEarth*t).*r(2,:);r(3,:)];
F.longitude_deg = @(x,y) atan2d(y,x);
F.geodetic_latitude_deg = @orbit_latitude;
F.noise_power_W = @(T_K,B_Hz) k*T_K.*B_Hz;
F.noise_power_dBW = @(T_K,B_Hz) 10*log10(k*T_K.*B_Hz);
F.GT_dBK = @(gain_dBi,T_K) gain_dBi-10*log10(T_K);
F.CN0_dBHz = @(CN_dB,B_Hz) CN_dB+10*log10(B_Hz);
F.CN_from_CN0_dB = @(CN0_dBHz,B_Hz) CN0_dBHz-10*log10(B_Hz);
F.CN0_link_dBHz = @(EIRP,loss,GT) EIRP-loss+GT-10*log10(k);
F.noise_temperature_K = @(noiseFigure_dB,T0_K) (10.^(noiseFigure_dB/10)-1).*T0_K;
F.system_temperature_K = @(antenna_K,receiver_K) antenna_K+receiver_K;
if nargout==0
    disp('Formula functions: positive linear inputs; gains and losses are linear unless named dB.');
    disp('Real detector SNR = A^2/noise variance; pulse matrices use rows for range and columns for pulses.');
    disp('Orbit inputs: km, radians, seconds; anomaly and time arrays are row vectors.');
    disp('Polar azimuth: counterclockwise from +X. Sampled Doppler interval: [-PRF/2, PRF/2).');
    disp(F);
end
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
