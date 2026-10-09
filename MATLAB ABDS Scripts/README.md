# MATLAB scripts for the 1090 MHz noise floor study

This is a reading copy of the 35 scripts, sorted into folders. The copy that runs lives in
`SDR Control/data`, where all 35 sit in one folder so MATLAB can find them.

## Folders

| Folder | Files | What's in it |
|---|---|---|
| (top level) | 1 | `expConfig.m`, the settings every script reads first |
| `1 Audit` | 2 | Checks our numbers against dump1090's, then compares its reported noise with the background |
| `2 Estimators` | 4 | Runs every noise estimator on real captures, made-up noise, and noise with messages added |
| `3 Detections` | 7 | Preamble counts against background, and how often messages decode at each noise level |
| `4 Figures` | 1 | Draws the four figures |
| `Helpers` | 17 | Functions the scripts call |
| `Tests` | 3 | Small checks of single functions |

## Run order

Each script needs tables made by the ones above it.

| # | Script | Folder | Tables it makes |
|---|---|---|---|
| 1 | `checkStep1` | 1 Audit | `step1.csv` |
| 2 | `auditMetrics` | 1 Audit | `audit.csv` |
| 3 | `estimatorsReal` | 2 Estimators | `real_est.csv`, `real_acf.csv` |
| 4 | `runSynth` | 2 Estimators | `synth_est.csv`, `synth_dump.csv` |
| 5 | `maskRecount` | 2 Estimators | `mask_noise.csv`, `mask_tail.csv` |
| 6 | `runOccupancy` | 2 Estimators | `occupancy.csv` |
| 7 | `fixCompare` | 3 Detections | `fix_compare.csv` |
| 8 | `runColorTwins` | 3 Detections | `color_twins.csv` |
| 9 | `runDither` | 3 Detections | `dither.csv` |
| 10 | `runPreambleLaw` | 3 Detections | `preamble_law.csv`, `preamble_law_fit.csv` |
| 11 | `runInjection` | 3 Detections | `injection.csv`, `injection_summary.csv` |
| 12 | `snrAnchor` | 3 Detections | `snr_anchor_cal.csv`, `snr_anchor_real.csv` |
| 13 | `makeFigures` | 4 Figures | the four figure files |

`falseAccept` (in `3 Detections`) and the three tests stand alone and can run any time.
