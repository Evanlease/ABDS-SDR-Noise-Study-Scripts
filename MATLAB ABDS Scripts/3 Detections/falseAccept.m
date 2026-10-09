% First we create a file of pure noise, the same one runSynth made at sigma
% 0.7174.

cfg = expConfig();
[bI, bQ] = synthNoise(0.7174, cfg.nSamples, 7174);
writeUC8(cfg.work, bI, bQ);

% Then we run dump1090 on it and start a for loop that prints every message it
% thought it found.

msg = runDump1090(cfg.work, 'raw');

for k = 1:numel(msg.hex)
    fprintf('%12d  %s  DF %d\n', msg.ts(k), msg.hex(k), bitshift(hex2dec(extractBefore(msg.hex(k), 3)), -3));
end
