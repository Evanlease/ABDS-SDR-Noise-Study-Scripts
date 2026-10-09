NOISE FLOOR OF AN 8-BIT ADS-B RECEIVER
MATLAB analysis code

Adrian Iugan
Captures taken September 21, 2026


READ THIS FIRST

This is the analysis code for a small study I ran on a 1090 MHz receiver. It's
organized for reading. It doesn't run from this layout, and the raw captures
aren't included. "Why it doesn't run from here" below says what's missing.

Everything comes from one receiver, at one site, on one day: seven captures of
30 seconds each. Read every number below as a measurement of my station, not as
a fact about RTL-SDR receivers in general.

All power levels are in dBFS (decibels relative to the largest signal the
receiver can record). They are uncalibrated. None of them are dBm.


THE QUESTION

ADS-B is the system where aircraft broadcast their own position on 1090 MHz.
Lots of people receive those broadcasts with a cheap USB radio (an RTL-SDR) and
a decoder program called dump1090-fa. The decoder reports a "noise" number, and
people graph it to judge how well their station is working.

I wanted to know three things:

  1. Does that noise number measure the receiver's background, or does
     aircraft traffic move it?
  2. Is there a way to measure the background that traffic doesn't move?
  3. What does a small background cost in decoded messages?


TERMS

  Background   The noise the receiver sees when no aircraft signal is there.

  LSB          One step of the receiver's 8-bit converter. It stores every
               sample as a whole number from 0 to 255, so one LSB is the
               smallest change it can record.

  Sigma        The spread (standard deviation) of the background, in LSB. At
               my station it's 0.71 to 0.89 LSB, which is less than one step.
               That's why rounding matters here.

  I and Q      The two numbers the receiver records for every sample.

  Mode S       The message format ADS-B uses. Mode A/C is an older, shorter
               reply that aircraft also send on 1090 MHz.

  Preamble     The short pulse pattern that starts every Mode S message. The
               decoder counts every spot where it thinks it sees one.

  SNR          Signal-to-noise ratio, in dB.

  Twin         A file of made-up noise built to match a real capture.


THE SETUP

  Receiver      RTL-SDR Blog V3 (RTL2832U chip with an 8-bit converter,
                R820T2 tuner)
  Antenna       The telescopic dipole that came with it, on a rooftop. I
                didn't record its length, height or orientation.
  Gain          Fixed at 44.5 dB
  Sample rate   2.4 million samples per second, 8-bit I and 8-bit Q
  Recorder      rtl_sdr
  Decoder       dump1090-fa 11.1, run on the recorded files
  Analysis      MATLAB R2026a, no toolboxes

The seven captures:

  #     Local time   Tuned to    Role      Sigma (LSB)   Mode S decoded
  001   09:25        1090 MHz    traffic   0.849         1,844
  002   13:13        1090 MHz    traffic   0.767         1,933
  003   16:27        1090 MHz    traffic   0.710         1,591
  004   16:57        1082 MHz    control   0.889         0
  005   22:02        1090 MHz    traffic   0.721         1,412
  006   22:03        1082 MHz    control   0.717         0
  007   22:04        1090 MHz    traffic   0.718         912

Each capture is 30 seconds, which is 72,000,000 samples and 144 MB.

Captures 004 and 006 are the controls. I tuned them to 1082 MHz, where a quick
scan beforehand found no loud pulses, so they show the background without
ADS-B traffic. Captures 005, 006 and 007 were taken within two minutes of each
other, so one control sits between two traffic captures.


WHAT I FOUND

1. The noise number dump1090 reports includes traffic.

   On the five traffic captures it sat 0.24 to 0.98 dB above the background.
   On the two controls it sat 0.00 and 0.06 dB above it. Across the three
   captures taken within two minutes, the reported number moved 1.0 dB while
   the background moved less than 0.05 dB.

2. Most of that extra power comes from signals the decoder never decoded.

   I hid every decoded message (Mode S and Mode A/C), plus 500 samples (208
   microseconds) on each side, and measured again. That removed only 14 to 44%
   of the extra power.

3. Two estimators that only look at the center of the converter's range held
   steady.

   The first counts the share of samples that land in the four center squares
   (the cell fraction). The second fits a rounded bell curve to the middle
   bytes (the central fit). On the real captures, both stayed within 0.04 dB
   of the control taken in the same two minutes. On made-up noise with no
   traffic, from 0.5 to 6 LSB, both were within 0.02 dB of the true value.

   They resist traffic. They aren't immune to it. In a made-up test at 0.72
   LSB, with test messages covering about 18% of the samples, the cell
   fraction read 1.02 dB high and the central fit read 0.46 dB high.
   dump1090's number rose 2.1 dB in the same test. My real captures had under
   1% of their samples covered by decoded messages.

4. The median-based estimators got stuck on the converter's steps.

   The clearest case is the statistic dump1090's automatic gain uses, a 40th
   percentile. It read -44.99 dBFS on all seven captures while the background
   power changed by 1.7 dB.

5. The preamble count follows the background.

   The count ran from 342,520 to 498,369 per capture. A straight line against
   sigma fit all seven captures within 1.8% (r = 0.998). Noise alone produces
   most of these: the two controls had no Mode S messages and still gave
   359,962 and 498,369.

6. A small background costs sensitivity. This one is a made-up test.

   I added test messages to made-up noise at 11 background levels, 900
   messages at each 1 dB step of SNR, and counted how many dump1090 decoded.
   At 0.72 LSB, which is about my station, messages needed 0.82 dB more SNR to
   decode half the time than they did at 6 LSB. At 0.5 LSB they needed 1.47 dB
   more.


WHAT I'M NOT CLAIMING

I'm not claiming any of this is new. Forum posts already pointed out that
dump1090's noise number includes traffic (obj in 2016, wiedehopf in 2020), and
the cell fraction idea comes from radio astronomy (Jenet and Anderson, 1998).
What this adds is the size of the effect at one station. I haven't finished my
own search of earlier work, so I can't say nobody has measured it before.

I'm not claiming the error is the same size anywhere else, or that any gain
setting is best.

Limits you should know about:

  - One dongle, one antenna, one site, one day, 3.5 minutes of data in total.
  - My site is quiet: 30 to 64 Mode S messages per second. What happens at a
    busy site is simulated only.
  - Every result above 0.89 LSB is made up, or made by adding noise to a real
    capture.
  - The 1082 MHz control isn't guaranteed to be empty. Other aviation systems
    use that band, and the edges of a strong 1090 MHz signal can reach it. I
    didn't test for either one.
  - Only dump1090-fa 11.1 was tested.


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

That's 35 MATLAB files. Folders 1 to 4 are numbered in the order they run.


RUN ORDER

Each script needs tables made by the ones above it. This order comes from
which tables each script reads and writes.

   1  checkStep1       (1 Audit)
      Four checks that the MATLAB math matches dump1090: its counts, the byte
      statistics, its noise number rebuilt from the raw samples, and the file
      reader and writer.
      Makes: step1.csv

   2  auditMetrics     (1 Audit)
      Finds each capture's background from the center squares and compares it
      with the noise dump1090 reported. Also works out what the automatic gain
      statistic would read.
      Makes: audit.csv

   3  estimatorsReal   (2 Estimators)
      Runs all eight estimators on each capture three ways: the whole file,
      the file with decoded messages cut out, and each 1 second block.
      Makes: real_est.csv, real_acf.csv

   4  runSynth         (2 Estimators)
      Runs the estimators on made-up noise at 21 sigmas from 0.5 to 6 LSB,
      where the true answer is known. Runs dump1090 on each file too.
      Makes: synth_est.csv, synth_dump.csv

   5  maskRecount      (2 Estimators)
      Hides decoded messages with wider and wider guards and measures how much
      extra power is left.
      Makes: mask_noise.csv, mask_tail.csv

   6  runOccupancy     (2 Estimators)
      Adds test messages to made-up noise until they fill 0 to 20% of the
      file, and tracks how far each estimator and dump1090's number rise.
      Makes: occupancy.csv

   7  fixCompare       (3 Detections)
      Runs dump1090 on each capture with and without --fix, its option for
      repairing messages with one wrong bit.
      Makes: fix_compare.csv

   8  runColorTwins    (3 Detections)
      Builds two twins of each capture and compares their preamble counts
      with the real one. One twin is plain noise. The other also copies how
      alike neighboring samples are in the real receiver.
      Makes: color_twins.csv

   9  runDither        (3 Detections)
      Adds noise to captures 003 and 006 and to their twins, to raise the
      background to 1.5, 2, 3 and 4 LSB.
      Makes: dither.csv

  10  runPreambleLaw   (3 Detections)
      Fits the straight line of preamble count against sigma, for the real
      captures and for their twins.
      Makes: preamble_law.csv, preamble_law_fit.csv

  11  runInjection     (3 Detections)
      Adds test messages at 23 SNR levels to noise at 11 backgrounds and
      records the share dump1090 decodes. Repeats it on control 006 and its
      twin as a check.
      Makes: injection.csv, injection_summary.csv

  12  snrAnchor        (3 Detections)
      Measures the SNR of every real Mode S message dump1090 decoded, after
      checking the method on made-up messages of known SNR.
      Makes: snr_anchor_cal.csv, snr_anchor_real.csv

  13  makeFigures      (4 Figures)
      Draws the four figures from the tables above.
      Makes: fig1_audit.png      reported noise next to the background
             fig2_preambles.png  preamble count against the background
             fig3_injection.png  how often messages decode against SNR
             fig4_occupancy.png  how far each estimate rises with traffic

falseAccept (3 Detections) and the three tests stand alone and can run any
time. falseAccept runs dump1090 on pure noise and prints every message it
thought it found.


THE EIGHT ESTIMATORS

All eight are in Helpers/allEstimators.m. Each one turns a set of samples into
a sigma.

  meanPower     The noise dump1090 reports, turned into a sigma
  sheppard      The spread of the bytes, minus the small amount rounding adds
  medianMag     From the middle value of the signal strength
  mad           From the typical distance of each byte from its middle value
  p40           From the 40th percentile that dump1090's automatic gain uses
  cellDC        The cell fraction, using the receiver's measured center
  cellCentred   The cell fraction, assuming the center is exactly 127.5
  centralFit    The central fit


THE HELPERS

  adaptiveP40       Copies the statistic dump1090's automatic gain uses
  allEstimators     Runs all eight estimators
  byteHist          Counts how many samples land on each pair of I and Q bytes
  cellSigma         Turns the share of samples in the center squares into sigma
  colorTwin         Builds a twin whose neighboring samples are as alike as
                    the real receiver's
  ditherIQ          Adds random noise to every byte, then rounds it
  dumpLUT           Rebuilds dump1090's table of signal strength for every
                    pair of bytes
  injectFrames      Adds test messages to a file of noise
  modesFrame        Builds one Mode S test message and its 24 check bits
  preambleLaw       Fits the straight line of preamble count against sigma
  quantModel        The chance of each byte when a bell curve gets rounded
  readDecoded       Reads the table of decoded messages for one capture
  readUC8           Reads a capture file into I and Q bytes
  runDump1090       Runs dump1090 on a file and pulls out its counts or its
                    messages
  sigmaFromCellDC   The cell fraction when the center isn't exactly 127.5
  synthNoise        Makes made-up noise at a chosen sigma
  writeUC8          Writes I and Q bytes back out as a capture file


THE TESTS

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
been run from this layout.