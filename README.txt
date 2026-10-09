NOISE FLOOR OF AN 8-BIT ADS-B RECEIVER
MATLAB analysis code

Adrian Iugan
Captures taken September 21, 2026


READ THIS FIRST

This folder holds the MATLAB code that did the analysis for a small study I
ran on a 1090 MHz receiver. This file is a guide to that code: how it's
organized, what each file does, and how the files depend on each other. It
only gives a short overview of the research.

The code is organized for reading. It doesn't run from this layout, and the
raw captures aren't included. "Why it doesn't run from here" at the bottom
says what's missing.


THE STUDY IN BRIEF

ADS-B is the system where aircraft broadcast their own position on 1090 MHz.
Lots of people receive it with a cheap USB radio (an RTL-SDR) and a decoder
program called dump1090-fa. The decoder reports a "noise" number that people
use to judge how well their station is working.

I recorded seven captures of 30 seconds each, on one receiver, at one site, on
one day. Five were on 1090 MHz with aircraft traffic. Two were controls tuned
to 1082 MHz, away from the ADS-B signal. The code asks three things of them:

  1. Does that noise number measure the receiver's background, or does
     aircraft traffic move it?
  2. Is there a way to measure the background that traffic doesn't move?
  3. What does a small background cost in decoded messages?

The short version: the reported number sat 0.24 to 0.98 dB above the
background on the traffic captures, and within 0.06 dB of it on the controls.
Two estimators that only look at the center of the converter's range stayed
within 0.04 dB of the control taken in the same two minutes. In a test with
made-up signals, messages at my station's background needed 0.82 dB more SNR
to decode half the time than they did at a large background.

These are measurements of one station. Power levels are uncalibrated dBFS, not
dBm. I'm not claiming any of it is new: forum posts had already pointed out
that the noise number includes traffic, and I haven't finished my own search
of earlier work.


WORDS USED BELOW

  Background   The noise the receiver sees when no aircraft signal is there.
  LSB          One step of the receiver's 8-bit converter, which stores every
               sample as a whole number from 0 to 255.
  Sigma        The spread (standard deviation) of the background, in LSB. At
               my station it's under one LSB, which is why rounding matters.
  Estimator    A formula that turns a set of samples into a sigma.
  I and Q      The two numbers the receiver records for every sample.
  Mode S       The message format ADS-B uses.
  Preamble     The short pulse pattern that starts every Mode S message.
  SNR          Signal-to-noise ratio, in dB.
  Twin         A file of made-up noise built to match a real capture.
  Table        A CSV file of results.


HOW THE CODE IS ORGANIZED

There are 35 MATLAB files and about 1,400 lines of code. They come in three
kinds.

  Settings (1 file)
      expConfig.m is a function that returns every setting in one place.

  Scripts (17 files)
      Each script runs one experiment from start to finish. The 14 in folders
      1 to 4 are the analysis. The 3 in Tests are small checks.

  Helpers (17 files)
      Each helper is a function that does one small job, like reading a
      capture file or making noise. 12 of the 17 are under 20 lines of code.

Four rules hold it together:

  1. Every analysis script starts with the same line, cfg = expConfig();
     which loads the settings.
  2. No script calls another script. Scripts only call helpers.
  3. Scripts pass results to each other through tables. A script reads the
     tables it needs, does its work, and saves its own. That's why the run
     order is fixed.
  4. Helpers never touch the result tables. Only scripts read and save those.

Most scripts have the same shape: load the settings, loop over the captures or
over a list of noise levels, collect the results in a table, then save the
table and print it. Each block of code has a short comment above it that says
what the block does.


HOW THE DATA MOVES

      The seven captures                Made-up noise and test
      and their tables of               messages, built by the
      decoded messages                  helpers
              |                                 |
              v                                 v
        The scripts in folders 1 to 3, which also run dump1090
                                |
                                v
                        18 result tables
                                |
                                v
                     makeFigures, in folder 4
                                |
                                v
                          four figures

Of the 14 analysis scripts, nine open the real captures. Three build made-up
noise instead (runSynth, runOccupancy and falseAccept). Two work only from
tables (runPreambleLaw and makeFigures). Nine of the 14 run dump1090.

A script that builds a made-up capture writes it to one scratch file and runs
dump1090 on that file, the same way it runs on a real capture.


WHAT'S IN EACH FOLDER

  (top level)     1 file     expConfig.m, the settings every script reads
                             first
  1 Audit         2 files    Checks my numbers against dump1090's, then
                             compares its reported noise with the background
  2 Estimators    4 files    Runs every noise estimator on real captures, on
                             made-up noise, and on noise with messages added
  3 Detections    7 files    Preamble counts against the background, and how
                             often messages decode at each noise level
  4 Figures       1 file     Draws the four figures
  Helpers         17 files   Functions the scripts call
  Tests           3 files    Small checks of single functions

Folders 1 to 4 are numbered in the order they run.


THE SETTINGS FILE

expConfig.m returns one variable, cfg, that holds:

  - Where things live: the captures, the tables of decoded messages, the
    result tables, and the scratch file for made-up captures
  - The path to dump1090 (line 15)
  - The sample rate (2.4 million per second), the samples in one capture (72
    million), and dump1090's clock rate (12 MHz)
  - The center and scale dump1090 uses to turn bytes into signal strength
  - The length of a preamble (19 samples)
  - The names of the seven captures, their short labels, and which five were
    on 1090 MHz


RUN ORDER

Each script needs tables made by the ones above it. "Reads" is what the script
takes in, "Calls" is the helpers it uses, and "Makes" is what it saves.

   1  checkStep1       (1 Audit)
      Four checks that the MATLAB math matches dump1090: its counts, the byte
      statistics, its noise number rebuilt from the raw samples, and the file
      reader and writer.
      Reads:  the captures and their decoded messages
      Calls:  byteHist, dumpLUT, readDecoded, readUC8, runDump1090, writeUC8
      Makes:  step1.csv

   2  auditMetrics     (1 Audit)
      Finds each capture's background from the center squares and compares it
      with the noise dump1090 reported. Also works out what the automatic gain
      statistic would read.
      Reads:  the captures and their decoded messages, step1.csv
      Calls:  adaptiveP40, byteHist, dumpLUT, quantModel, readDecoded,
              readUC8, sigmaFromCellDC
      Makes:  audit.csv

   3  estimatorsReal   (2 Estimators)
      Runs all eight estimators on each capture three ways: the whole file,
      the file with decoded messages cut out, and each 1 second block.
      Reads:  the captures and their decoded messages, audit.csv
      Calls:  allEstimators, byteHist, readDecoded, readUC8, synthNoise
      Makes:  real_est.csv, real_acf.csv

   4  runSynth         (2 Estimators)
      Runs the estimators on made-up noise at 21 sigmas from 0.5 to 6 LSB,
      where the true answer is known. Runs dump1090 on each file too.
      Reads:  audit.csv
      Calls:  allEstimators, byteHist, dumpLUT, quantModel, runDump1090,
              synthNoise, writeUC8
      Makes:  synth_est.csv, synth_dump.csv

   5  maskRecount      (2 Estimators)
      Hides decoded messages with wider and wider guards and measures how much
      extra power is left.
      Reads:  the captures and their decoded messages, audit.csv, real_acf.csv
      Calls:  colorTwin, dumpLUT, quantModel, readDecoded, readUC8
      Makes:  mask_noise.csv, mask_tail.csv

   6  runOccupancy     (2 Estimators)
      Adds test messages to made-up noise until they fill 0 to 20% of the
      file, and tracks how far each estimator and dump1090's number rise.
      Reads:  real_est.csv, real_acf.csv
      Calls:  allEstimators, byteHist, colorTwin, dumpLUT, injectFrames,
              modesFrame, quantModel, runDump1090, writeUC8
      Makes:  occupancy.csv

   7  fixCompare       (3 Detections)
      Runs dump1090 on each capture with and without --fix, its option for
      repairing messages with one wrong bit.
      Reads:  the captures
      Calls:  runDump1090
      Makes:  fix_compare.csv

   8  runColorTwins    (3 Detections)
      Builds a twin of each capture that also copies how alike neighboring
      samples are in the real receiver. Compares its preamble count with the
      real capture's and with the plain noise from runSynth.
      Reads:  the captures, audit.csv, real_est.csv, real_acf.csv,
              synth_dump.csv
      Calls:  colorTwin, runDump1090, writeUC8
      Makes:  color_twins.csv

   9  runDither        (3 Detections)
      Adds noise to captures 003 and 006 and to their twins, to raise the
      background to 1.5, 2, 3 and 4 LSB.
      Reads:  the captures, real_est.csv, real_acf.csv
      Calls:  allEstimators, byteHist, colorTwin, ditherIQ, readUC8,
              runDump1090, writeUC8
      Makes:  dither.csv

  10  runPreambleLaw   (3 Detections)
      Fits the straight line of preamble count against sigma, for the real
      captures and for their twins.
      Reads:  audit.csv, fix_compare.csv, color_twins.csv, dither.csv
      Calls:  preambleLaw
      Makes:  preamble_law.csv, preamble_law_fit.csv

  11  runInjection     (3 Detections)
      Adds test messages at 23 SNR levels to noise at 11 backgrounds and
      records the share dump1090 decodes. Repeats it on control 006 and its
      twin as a check.
      Reads:  the captures, real_est.csv, real_acf.csv
      Calls:  colorTwin, injectFrames, modesFrame, readUC8, runDump1090,
              writeUC8
      Makes:  injection.csv, injection_summary.csv

  12  snrAnchor        (3 Detections)
      Measures the SNR of every real Mode S message dump1090 decoded, after
      checking the method on made-up messages of known SNR.
      Reads:  the captures and their decoded messages, real_est.csv,
              real_acf.csv, injection.csv
      Calls:  colorTwin, injectFrames, modesFrame, readDecoded, readUC8,
              runDump1090, writeUC8
      Makes:  snr_anchor_cal.csv, snr_anchor_real.csv

  13  makeFigures      (4 Figures)
      Draws the four figures from 7 of the 18 tables.
      Reads:  audit.csv, fix_compare.csv, color_twins.csv, dither.csv,
              injection.csv, injection_summary.csv, occupancy.csv
      Calls:  dumpLUT, preambleLaw, quantModel
      Makes:  fig1_audit.png      reported noise next to the background
              fig2_preambles.png  preamble count against the background
              fig3_injection.png  how often messages decode against SNR
              fig4_occupancy.png  how far each estimate rises with traffic

falseAccept (3 Detections) stands alone. It runs dump1090 on pure noise and
prints every message it thought it found. It calls runDump1090, synthNoise and
writeUC8, and it saves nothing.


THE HELPERS

I've grouped them here by the job they do. In the folder they're one flat list.

  Reading and writing files

    readUC8           Reads a capture file into I and Q bytes
    writeUC8          Writes I and Q bytes back out as a capture file
    readDecoded       Reads the table of decoded messages for one capture

  Working with dump1090

    runDump1090       Runs dump1090 on a file and pulls out its counts or its
                      messages
    dumpLUT           Rebuilds dump1090's table of signal strength for every
                      pair of bytes
    adaptiveP40       Copies the statistic dump1090's automatic gain uses

  Estimating the noise

    allEstimators     Runs all eight estimators
    byteHist          Counts how many samples land on each pair of I and Q
                      bytes
    cellSigma         Turns the share of samples in the center squares into
                      sigma
    sigmaFromCellDC   The same, when the center isn't exactly 127.5
    quantModel        The chance of each byte when a bell curve gets rounded

  Making test signals

    synthNoise        Makes made-up noise at a chosen sigma
    colorTwin         Builds a twin whose neighboring samples are as alike as
                      the real receiver's
    ditherIQ          Adds random noise to every byte, then rounds it
    modesFrame        Builds one Mode S test message and its 24 check bits
    injectFrames      Adds test messages to a file of noise

  Fitting

    preambleLaw       Fits the straight line of preamble count against sigma

Helpers mostly stand alone. allEstimators is the exception: it calls cellSigma,
dumpLUT and sigmaFromCellDC. Four helpers load the settings themselves
(adaptiveP40, dumpLUT, readDecoded and runDump1090).


THE EIGHT ESTIMATORS

All eight are in Helpers/allEstimators.m. Each one turns a set of samples into
a sigma.

  meanPower     The noise dump1090 reports, turned into a sigma
  sheppard      The spread of the bytes, minus the small amount rounding adds
  medianMag     From the middle value of the signal strength
  mad           From the typical distance of each byte from its middle value
  p40           From the 40th percentile that dump1090's automatic gain uses
  cellDC        The share of samples in the four center squares, using the
                receiver's measured center
  cellCentred   The same, assuming the center is exactly 127.5
  centralFit    Fits a rounded bell curve to the middle bytes


THE TESTS

Each test checks one helper against numbers typed into the test. They don't
need the captures or dump1090.

  testCellSigma     Checks cellSigma against the two controls
  testModesFrame    Checks modesFrame against a published example message
  testPreambleLaw   Checks preambleLaw's slope, intercept and r


WHAT YOU'D NEED TO RUN IT

  - MATLAB. I used R2026a with no toolboxes. I haven't tried older versions.
  - dump1090-fa. I used version 11.1 from Homebrew on a Mac. Line 15 of
    expConfig.m holds the path to it, so change that line for your machine.
  - The seven capture files. They're raw 8-bit I and Q, the format rtl_sdr
    writes.
  - One table of decoded messages per capture, with the sample where each
    message starts and ends. A separate script made those, and it isn't in
    this folder.

The made-up noise uses fixed random seeds, so those files come out the same
every time.


WHY IT DOESN'T RUN FROM HERE

Two reasons.

First, MATLAB only finds a function when its file is in the current folder or
on its path (the list of folders MATLAB searches). I sorted the files into
folders so they're easier to read, and nothing adds those folders to the path.
Putting all 35 files back in one folder fixes this.

Second, expConfig.m looks for the captures in a folder called iq next to
itself, writes results to a folder called out next to itself, and looks for
the decoded tables in a folder one level up. None of those are here.

Everything was run from one flat folder with the data beside it. It has not
been run from this layout. The run order above comes from which tables each
script reads and writes.
