% First we start a for loop that runs dump1090 on each capture twice, without
% --fix and with it.

cfg = expConfig();
R = table();

for c = 1:numel(cfg.names)
    f = fullfile(cfg.iqDir, [cfg.names{c} '.iq']);
    a = runDump1090(f, 'stats');  b = runDump1090(f, 'stats', '--fix');

    R = [R; table(string(cfg.labels{c}), a.accepted, b.accepted, b.repaired, a.modeac, b.modeac, ...
                  a.noiseDbfs, b.noiseDbfs, a.preambles, b.preambles, ...
                  'VariableNames', {'capture', 'accepted', 'acceptedFix', 'repairedFix', 'modeac', 'modeacFix', ...
                                    'noise', 'noiseFix', 'preambles', 'preamblesFix'})];

end

% Then we find the percent of extra messages that --fix gets us.

R.gainPct = 100 * ((R.acceptedFix + R.repairedFix) ./ R.accepted - 1);
R.gainPct(R.accepted == 0) = NaN;

% Lastly, we save and print the table.

writetable(R, fullfile(cfg.outDir, 'fix_compare.csv'));
disp(R)
