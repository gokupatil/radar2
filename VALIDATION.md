# Review and validation

Prepared from the uploaded 23-page `Radar all codes.pdf` and the nine-experiment index image. Source parameters were retained where meaningful. The files are corrected educational simulations, with figures arranged for readability rather than exact replicas of the PDF screenshots.

## What was tested

- All 11 `.m` files parsed without grammar errors or missing syntax nodes using tree-sitter-matlab 1.3.4.
- The syntax trees contain zero comment nodes. No explanatory comments or section comments occur in the MATLAB files.
- Independent numerical implementations in Python with NumPy and SciPy passed the checks below.
- The combined file embeds the same experiment implementations as the individual files, including the orbit helpers.
- The archive contains exactly 11 `.m` files plus `README.md` and this validation record.

**Not tested:** execution in MATLAB or Octave, MATLAB's Code Analyzer, actual MATLAB plotting/animation, GitHub upload, or cloning into your MATLAB installation. Neither MATLAB nor Octave was available here. A successful independent calculation is not an actual execution test of a `.m` file. Python and MATLAB random generators need not produce identical samples; Monte Carlo values should agree with theory within statistical variation.

## Per-experiment review

| No. | Source review and change | Independent numerical result / expected MATLAB behavior |
| --- | --- | --- |
| 1 | The source radar equation is correct. Retained transmit power, gains, RCS, losses and sensitivity; used explicit vector operations, labeled units and added checks. The frequency sweep holds gain constant. | Maximum range **10.3312585 km**. Received power at that range **-100 dBm**. Doubling range reduces received power by 16. Five curves show power decay and the expected parameter trends. |
| 2 | The source `awgn(...,'measured')` calls use a unit signal under H1 and a 0.1 baseline under H0, producing different noise powers and a nonzero noise-only mean. Replaced with zero-mean H0 and identical noise variance under both hypotheses. Removed the toolbox dependency. Swept Pfa to generate full useful ROC curves and compared Monte Carlo with exact real-Gaussian theory. | At 0 dB and threshold 0.5: **Pd = 0.69146246**, **Pfa = 0.30853754**. Independent 200,000-trial estimates differed from theory by at most **0.002856** across seven SNRs. ROC probabilities increase with SNR and Pfa. The original false-alarm screenshot should change. |
| 3 | Directly sampling a 10 GHz carrier at 100 kHz aliases the carrier. Replaced it with the equivalent complex baseband Doppler signal; retained the velocity and carrier-frequency sweeps. Corrected PSD normalization and used a signed frequency axis. | A 60 m/s target at 10 GHz gives **4000 Hz**, recovered by a 100 Hz-bin FFT. No attempt to sample the RF carrier is needed. |
| 4 | The original zero-padded delay creates a startup transient. Used adjacent valid pulse differences and corresponding time alignment. Added analytical frequency response. | Constant clutter cancels to **zero**. At 100 Hz Doppler and 1000 Hz PRF, amplitude gain is **0.61803399**. The canceller attenuates the moving target too; nulls occur at integer multiples of PRF. |
| 5 | Used exactly 20 pulse samples instead of an endpoint-inclusive 21. Replaced circular correlation with zero-padded linear matched filtering and corrected its delay indexing. Replaced endpoint-inclusive Doppler `linspace` with exact FFT bins. Used `20*log10` for amplitude, a Doppler window, and an upward range axis. Added a fast-time Doppler phase and checks on estimated target locations. | Peaks at **3000 m, 11.953125 m/s** and **6000 m, -7.968750 m/s**, within one velocity bin of 12 and -8 m/s. Resolution **30 m**; sampled range spacing **7.5 m**; velocity-bin spacing **0.234375 m/s**. Direct convolution and FFT matched filtering agreed within 1e-12. |
| 6 | Reconstructed the display initialization and update logic from the scanned source. Added reflected boundary handling, circular beam-angle differences and safe handling when the animation window is closed. A-scan now shows echoes only within the current beam. B-scan and red PPI markers are explicitly simulated target truth; green PPI markers indicate recent beam detections. | All **900** motion frames stayed within 0–100 km. Final ranges approximately **[15, 3, 86] km**. Every target was illuminated; independent detection-frame counts **[10, 10, 8]**. Cartesian coordinate norms equal target ranges. Graphics themselves require MATLAB verification. |
| 7 | Source link-budget arithmetic is correct. Preserved both link inputs, including **40 dBi** downlink receive gain. Added explicit link-closure status and combined C/N under independent-noise assumptions. | Uplink/downlink received power **-120.488311 / -122.159675 dBW**; C/N **4.511689 / 2.840325 dB**; margin **-5.488311 / -7.159675 dB**. These negative margins match the source and mean the supplied links fail the 10 dB requirement. Changing the arithmetic to force positive margins would be incorrect. |
| 8 | The supplied code calls an absent `keplerOrbitECI` helper and uses undefined longitude/latitude arrays. Added a Newton solver, orbital rotations, rotating-Earth ECI-to-ECEF conversion and iterative ellipsoidal geodetic latitude. Preserved the source orbital elements and antimeridian plot breaks. | Periods: LEO **1.576938 h**, MEO **11.978303 h**, GEO **23.934442 h**. Kepler residuals below 1e-11; radii within perigee/apogee limits. GEO latitude is zero and longitude varies by about **0.02269 degrees** over 24 h. |
| 9 | Source thermal-noise, C/N and G/T relationships are correct. Retained base inputs; used the exact Boltzmann constant **1.380649e-23 J/K** and added C/N0. | Noise **2.0709735e-15 W = -146.838255 dBW**; C/N **46.838255 dB**; G/T **13.239087 dB/K**; C/N0 **106.838255 dB-Hz**. Doubling temperature raises noise by **3.010300 dB**. |

## Model conventions and limits

- Experiment 2 is the source's **real-valued Gaussian amplitude detector**, with SNR defined as A^2 divided by real noise variance. It is not the default complex coherent model of `rocsnr`; those conventions differ. Rare false alarms are represented by the theoretical ROC, not claimed measured by 200,000 Monte Carlo trials.
- Positive velocity corresponds to positive Doppler in the signal model. Experiment 5 uses complex slow-time samples and unambiguous velocity interval [-15, 15) m/s. Targets use the usual fixed-range, short coherent-processing-interval approximation; it is a matched-filter demonstration rather than a CFAR detector. The peak search uses small neighborhoods around the simulated ranges as a validation check, not an unknown-target search algorithm.
- Experiment 6 uses mathematical azimuth (counterclockwise from +X), matching the source Cartesian conversion, rather than bearing clockwise from north. Motion rates are per animation frame, not a claimed physical radar revisit time. The display is a simulation, not measured radar returns or an estimated tracking filter.
- Experiment 7 takes a fixed noise power of -125 dBW for each link as supplied. It does not infer a temperature or bandwidth. Combined C/N assumes independent uplink and downlink noise referred to a common carrier and bandwidth, with no extra transponder/interference contribution.
- Experiment 8 is a two-body Keplerian model with Earth rotation. No drag, J2 perturbation, real epoch, GMST calculation or map dataset is used. The source GEO orbit has a small nonzero eccentricity and nominal 35786 km altitude, so slight drift is expected. A sidereal day is about 23.934 h, not 24 h. Ground-track latitude is the geodetic latitude of the satellite's ECEF position, using the WGS84 ellipsoid; it differs from spherical geocentric latitude.
- The formulas file includes reusable equations for all nine topics plus closely related noise-temperature and C/N0 relations. Physical input domains and array conventions are described in `README.md`.

## Finish validation in MATLAB

First run:

```matlab
results = all_experiments(false,false);
```

Check that no assertions fail and the completion message appears. Then run:

```matlab
results = all_experiments(true,false);
```

Compare the plots and printed values with the table. Finally check the two animations:

```matlab
exp06_radar_display(true,true);
exp08_satellite_orbits(true,true);
```

The MATLAB runtime checks cover key numerical invariants. Visual layout, refresh behavior and animation performance remain dependent on your MATLAB installation.
