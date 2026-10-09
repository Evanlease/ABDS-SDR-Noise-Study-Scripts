% First we load the control autocorrelation and each capture's centralFit sigma.

cfg = expConfig();
rest = readtable(fullfile(cfg.outDir, 'real_est.csv'));
acf = readtable(fullfile(cfg.outDir, 'real_acf.csv'));

ctl = ismember(string(acf.capture), ["004 control", "006 control"]);
target = mean([acf.II1(ctl), acf.II2(ctl), acf.II3(ctl)], 1);
sigCF = rest.estFull(string(rest.estimator) == "centralFit");

% Then we start a for loop over captures 3 and 6, where we add noise to the real
% file and its twin, run the estimators and dump1090, and count the loud samples.

T = table();

for c = [3 6]
    [rI, rQ] = readUC8(fullfile(cfg.iqDir, [cfg.names{c} '.iq']), 1, cfg.nSamples);
    [tI, tQ] = colorTwin(sigCF(c), cfg.nSamples, round(sigCF(c) * 1e4), target);

    for eff = [0 1.5 2 3 4]
        if eff == 0
            sd = 0;  effS = sqrt(sigCF(c)^2 + 1/12);
        else
            sd = sqrt(eff^2 - sigCF(c)^2 - 1/12);  effS = eff;
        end

        for src = ["real", "twin"]
            if src == "real", bI = rI; bQ = rQ; else, bI = tI; bQ = tQ; end

            if sd > 0
                [bI, bQ] = ditherIQ(bI, bQ, sd, round(eff * 100) + 7 * c);
            end

            [e, names] = allEstimators(byteHist(bI, bQ));
            writeUC8(cfg.work, bI, bQ);
            st = runDump1090(cfg.work, 'stats');

            z2 = (double(bI) - mean(double(bI))).^2 + (double(bQ) - mean(double(bQ))).^2;
            tails = arrayfun(@(k) sum(z2 > (k * effS)^2), [3 4 5]);
            ex = exp(-[3 4 5].^2 / 2) * numel(bI);

            T = [T; table(string(cfg.labels{c}), src, eff, sd, effS, ...
                          e(strcmp(names, 'cellDC')), e(strcmp(names, 'centralFit')), ...
                          e(strcmp(names, 'meanPower')), e(strcmp(names, 'sheppard')), ...
                          st.preambles, st.accepted, st.modeac, st.noiseDbfs, ...
                          tails(1) / ex(1), tails(2) / ex(2), tails(3), ex(3), ...
                          'VariableNames', {'capture', 'src', 'target', 'sigmaD', 'effSigma', ...
                                            'cellDC', 'centralFit', 'meanPower', 'sheppard', ...
                                            'preambles', 'accepted', 'modeac', 'noise', ...
                                            'tail3Ratio', 'tail4Ratio', 'tail5', 'tail5Rayleigh'})];

            fprintf('%s %s eff %.2f done\n', cfg.labels{c}, src, effS);

        end
    end
end

% Lastly, we save the table and print each estimator's error in dB.

writetable(T, fullfile(cfg.outDir, 'dither.csv'));

T.cellDCdb = 20 * log10(T.cellDC ./ T.effSigma);
T.centralFitDb = 20 * log10(T.centralFit ./ T.effSigma);
T.meanPowerDb = 20 * log10(T.meanPower ./ T.effSigma);

disp(T(:, {'capture', 'src', 'effSigma', 'cellDCdb', 'centralFitDb', 'meanPowerDb', ...
           'preambles', 'accepted', 'modeac', 'noise', 'tail3Ratio', 'tail4Ratio', 'tail5', 'tail5Rayleigh'}))
