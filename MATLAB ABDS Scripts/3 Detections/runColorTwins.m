% First we load earlier results and average the autocorrelation of the two
% controls.

cfg = expConfig();
aud = readtable(fullfile(cfg.outDir, 'audit.csv'));
acf = readtable(fullfile(cfg.outDir, 'real_acf.csv'));
syn = readtable(fullfile(cfg.outDir, 'synth_dump.csv'));

ctl = ismember(string(acf.capture), ["004 control", "006 control"]);
target = mean([acf.II1(ctl), acf.II2(ctl), acf.II3(ctl)], 1);

fprintf('target byte ACF lags 1-3: %.4f %.4f %.4f\n', target);

% Then we get each capture's centralFit sigma.

rest = readtable(fullfile(cfg.outDir, 'real_est.csv'));
sigCF = rest.estFull(string(rest.estimator) == "centralFit");

% Start a for loop over the captures, where we build two twins, run dump1090 on
% them, and compare the preamble counts with the real capture.

C = table();

for c = 1:numel(cfg.names)
    s = round(aud.sigmaDC(c), 4);
    real = runDump1090(fullfile(cfg.iqDir, [cfg.names{c} '.iq']), 'stats');

    [bI, bQ] = colorTwin(sigCF(c), cfg.nSamples, round(s * 1e4), target);
    writeUC8(cfg.work, bI, bQ);
    stCF = runDump1090(cfg.work, 'stats');

    [bI, bQ, h] = colorTwin(s, cfg.nSamples, round(s * 1e4), target);
    x = double(bI) - mean(double(bI));
    got = arrayfun(@(L) mean(x(1:end-L) .* x(1+L:end)), 1:3) / mean(x.^2);

    writeUC8(cfg.work, bI, bQ);
    st = runDump1090(cfg.work, 'stats');

    white = syn.preambles(abs(syn.sigma - s) < 5e-5);
    C = [C; table(string(cfg.labels{c}), s, real.preambles, white, st.preambles, ...
                  100 * (white / real.preambles - 1), 100 * (st.preambles / real.preambles - 1), ...
                  sigCF(c), 100 * (stCF.preambles / real.preambles - 1), ...
                  st.modeac, st.accepted, st.noiseDbfs, got(1), got(2), got(3), ...
                  'VariableNames', {'capture', 'sigma', 'realPre', 'whitePre', 'colorPre', ...
                                    'whiteExcessPct', 'colorExcessPct', 'sigmaCF', 'colorCFExcessPct', 'colorModeac', ...
                                    'colorAccepted', 'colorNoise', 'acf1', 'acf2', 'acf3'})];

end

% Lastly, we print the filter and save the table.

fprintf('FIR taps: %s\n', mat2str(h, 4));
writetable(C, fullfile(cfg.outDir, 'color_twins.csv'));
disp(C)
