function cfg = expConfig()

% First we set up where the code, the captures and the results live.

here = fileparts(mfilename('fullpath'));
cfg.root     = fileparts(here);
cfg.expDir   = here;

cfg.iqDir    = fullfile(here, 'iq');
cfg.decDir   = fullfile(cfg.root, 'paper-a-noise', 'decoded');
cfg.work     = fullfile(cfg.iqDir, 'derived', 'work.iq');

cfg.outDir   = fullfile(here, 'out');
cfg.dumpLog  = fullfile(cfg.outDir, 'dump1090_last.txt');
cfg.dump1090 = '/opt/homebrew/bin/dump1090';

% Then we set the sample rate, dump1090's clock rate, and the samples in a 30 s
% capture.

cfg.fs       = 2.4e6;
cfg.mlatHz   = 12e6;
cfg.nSamples = 72e6;

% These are the center and scale dump1090 uses to turn bytes into signal
% strength.

cfg.dcDump   = 127.4;
cfg.fsDump   = 128;

% The preamble is 19 samples long, so the data starts 19 samples in.

cfg.preambleSamples = 19;

% Next we list the seven captures, their short labels, and which were on 1090
% MHz.

cfg.names = {'20260921T132539Z_morning', '20260921T171320Z_midday', ...
             '20260921T202723Z_afternoon', '20260921T205711Z_control', ...
             '20260922T020220Z_late', '20260922T020320Z_control', ...
             '20260922T020354Z_late2'};

cfg.labels = {'001 morning', '002 midday', '003 afternoon', '004 control', ...
              '005 late', '006 control', '007 late2'};
cfg.is1090 = logical([1 1 1 0 1 0 1]);

% Lastly, we make the output and scratch folders if they don't exist yet.

if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
if ~exist(fileparts(cfg.work), 'dir'), mkdir(fileparts(cfg.work)); end

end
