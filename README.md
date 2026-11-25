# Radar and satellite communication MATLAB experiments

This package contains exactly 11 `.m` files: nine independent experiments, one self-contained combined implementation, and one reusable formula collection. No MATLAB file contains comments. Format tokens inside printed strings, such as `%.2f`, are not comments.

The implementations use base MATLAB functions and do not require Communications, Phased Array System, Radar, Aerospace, Mapping, or Satellite Communications Toolbox. MATLAB R2023b or newer is recommended for the Git cloning instructions. The numerical code does not depend on `gitclone`.

## Verification status

All 11 files passed a MATLAB grammar parser, including a check for absence of comments. Independent NumPy/SciPy calculations passed checks for all nine experiments. MATLAB and Octave were unavailable in the preparation environment, so the actual `.m` files, MATLAB graphics, and animations have **not** been executed here. See `VALIDATION.md` for corrections, expected outputs, and limits. The files contain runtime assertions for important invariants; these run when you execute them in MATLAB.

## Files

| File | Experiment / purpose |
| --- | --- |
| `exp01_radar_range.m` | Radar range equation and five parameter sweeps |
| `exp02_radar_detection_roc.m` | Detection in real Gaussian noise, Monte Carlo comparison, ROC curves |
| `exp03_doppler_shift.m` | Doppler frequency vs velocity/frequency and baseband FFT |
| `exp04_mti_canceller.m` | Two-pulse delay-line cancellation, spectra, frequency response |
| `exp05_pulse_doppler.m` | Two moving targets, matched filtering, range-Doppler map |
| `exp06_radar_display.m` | A-scan, B-scan, moving-target PPI animation |
| `exp07_satellite_link_budget.m` | Uplink/downlink powers, C/N, margins, distance sweep |
| `exp08_satellite_orbits.m` | LEO/MEO/GEO orbits, animation, 24-hour ground tracks |
| `exp09_system_noise.m` | Noise power, C/N, C/N0, G/T and four sweeps |
| `all_experiments.m` | All nine implementations and required local helpers in one file |
| `all_formulas.m` | Callable equations, constants and orbital helper functions |

## Run in MATLAB

Open the extracted folder as MATLAB's Current Folder. These are function files: call their names in the Command Window or open a file and click Run. Do not paste function definitions into the Command Window.

Run every experiment with figures and animations:

```matlab
results = all_experiments;
```

Run with figures but without animation:

```matlab
results = all_experiments(true,false);
```

Run numerical calculations and their internal checks without figures:

```matlab
results = all_experiments(false,false);
```

Run one experiment:

```matlab
r1 = exp01_radar_range;
r5 = exp05_pulse_doppler;
r6 = exp06_radar_display(true,false);
r8 = exp08_satellite_orbits(true,false);
```

All individual experiments accept a first argument `makePlots`. Experiments 6 and 8 also accept `animate`. Both flags default to `true`. `all_experiments` returns a 1-by-9 cell array of result structures, while each individual file returns its own structure. `all_experiments.m` embeds the complete implementations; it can run without the nine individual files. Random experiments use fixed seeds and restore your previous random state on exit. Existing figures are not closed automatically.

Example formula calculations:

```matlab
F = all_formulas;
lambda = F.wavelength_m(10e9)
fd = F.doppler_Hz(60,lambda)
noise_W = F.noise_power_W(150,1e6)
GT_dBK = F.GT_dBK(35,150)
```

Expected results: `0.03`, `4000`, `2.0709735e-15`, and approximately `13.239087`.

Use positive frequency, range, temperature, bandwidth and linear power inputs. Radar antenna gains and losses are linear ratios unless a parameter name explicitly says dB. Orbit helpers use km, seconds and radians, with time/anomaly samples in row vectors. `two_pulse_canceller` operates along columns (slow time). ECI-to-ECEF conversion assumes a zero initial Earth rotation angle. The formula collection includes the main physical equations and processing transforms used in the experiments, not UI/plotting operations.

## Upload to GitHub

1. Download and extract `radar_matlab.zip`.
2. Sign in to GitHub and open https://github.com/gokupatil/radar . If it does not exist, create a repository named `radar` under `gokupatil`.
3. Select **Add file > Upload files**. For a new empty repository, use **uploading an existing file**.
4. Upload the 11 `.m` files and these two Markdown documents from the extracted folder. Upload the files themselves into the repository root; uploading only the ZIP does not make the MATLAB functions available after cloning.
5. Enter a commit message such as `Add corrected radar and satellite experiments`.
6. Commit to the main branch if permitted. If you choose a new branch, create and merge its pull request before cloning the default branch.
7. Verify the repository's main page lists `all_experiments.m` and the other files.

## Clone into MATLAB

`gitclone` is available starting in MATLAB R2023b. Choose a writable Current Folder with no existing `radar` subfolder, then run:

```matlab
repo = gitclone("https://github.com/gokupatil/radar.git");
cd(repo.WorkingFolder);
results = all_experiments(true,false);
```

Your original command also works for the cloning step:

```matlab
gitclone("https://github.com/gokupatil/radar.git")
```

Saving its return value as `repo` makes changing to the cloned folder reliable. Cloning downloads the repository; the final command runs the experiments. The examples assume you uploaded the files at the root. If they are in a subfolder, navigate into it before running them.

For a private repository, authenticate with GitHub when MATLAB requests it. If a `radar` clone already exists, open that clone instead of cloning over it, or clone to a different empty folder. If `gitclone` is undefined, use MATLAB's Git clone interface or download **Code > Download ZIP** from GitHub and open the extracted folder in MATLAB.

After cloning, a quick numerical check is:

```matlab
results = all_experiments(false,false);
```

Success prints `All nine experiments completed and internal checks passed.` For full visual verification, run with figures and animations and compare with the expected behavior in `VALIDATION.md`.

References:
- https://www.mathworks.com/help/matlab/ref/gitclone.html
- https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository
- https://www.mathworks.com/help/phased/ref/rocsnr.html
- https://www.mathworks.com/help/satcom/gs/satellite-link-budget.html
