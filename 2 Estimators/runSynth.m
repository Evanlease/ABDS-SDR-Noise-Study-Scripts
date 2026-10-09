% First we make the list of sigmas to test, from 0.5 to 6 plus the seven
% measured ones.

cfg = expConfig();
aud = readtable(fullfile(cfg.outDir, 'audit.csv'));

sigmas = unique([0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 2 2.5 3 4 5 6, ...
                 round(aud.sigmaDC', 4)]);
nBlk = 30;  blkLen = 2.4e6;

% Then we create two empty tables and two checks that start out passing.

E = table();  D = table();
gateCell = true;  gateP40 = true;

% Start a for loop over the sigmas, where we make a 30 s file of made-up noise
% and run every estimator on it.

for s = sigmas
    seed = round(s * 1e4);
    [bI, bQ] = synthNoise(s, cfg.nSamples, seed);

    H = byteHist(bI, bQ);
    N = sum(H(:));
    [full, names] = allEstimators(H);

    % Then we run the estimators on each 1 s block too.

    blk = nan(nBlk, numel(names));

    for k = 1:nBlk
        idx = (k - 1) * blkLen + (1:blkLen);
        blk(k, :) = allEstimators(byteHist(bI(idx), bQ(idx)));
    end

    % Next we add a row for each estimator with its error in dB.

    for e = 1:numel(names)
        row = table("base", s, string(names{e}), full(e), ...
                    20 * log10(full(e) / s), mean(20 * log10(blk(:, e) / s)), ...
                    std(20 * log10(blk(:, e) / s)), ...
                    'VariableNames', {'variant', 'sigma', 'estimator', 'estFull', ...
                                      'biasDb', 'blockBiasDb', 'blockStdDb'});
        E = [E; row];
    end

    % After that we check the share of samples in the 4 center squares matches
    % the model.

    p = quantModel(127.397, s);
    cellModel = sum(p(128:129))^2;
    cellMeas = sum(sum(H(128:129, 128:129))) / N;

    se = sqrt(cellModel * (1 - cellModel) / N);
    if abs(cellMeas - cellModel) > 4 * se, gateCell = false; end

    % We also check the 40th percentile level matches the model.

    J = p * p';
    lut = dumpLUT();
    [ml, ~, g] = unique(lut(:));
    pm = accumarray(g, J(:));

    levModel = ml(find(cumsum(pm) > 0.40, 1));
    levMeas = full(5) * sqrt(-2 * log(0.6)) / 128 * 65535;
    if abs(levMeas - levModel) > 0.5, gateP40 = false; end

    % Finish off the loop by running dump1090 on the file and saving its counts.

    writeUC8(cfg.work, bI, bQ);
    st = runDump1090(cfg.work, 'stats');

    D = [D; table(s, st.preambles, st.badFormat, st.accepted, st.modeac, ...
                  st.noiseDbfs, 10 * log10(sum(H(:) .* reshape(dumpLUT(), [], 1).^2) / N / 65535^2), ...
                  100 * cellMeas, 100 * cellModel, ...
                  'VariableNames', {'sigma', 'preambles', 'bad', 'accepted', 'modeac', ...
                                    'noiseReported', 'totalDb', 'cellPct', 'cellModelPct'})];

    fprintf('sigma %.4f  preambles %d  noise %.1f  cell %.4f / %.4f\n', ...
            s, st.preambles, st.noiseDbfs, 100 * cellMeas, 100 * cellModel);

end

% Next we start a for loop that runs the estimators on two variations: centered
% noise, and Q 5% bigger than I.

for variant = ["centred", "imbalance"]
    for s = [0.6 0.8 1.0 1.5 2 4 6]
        if variant == "centred"
            [bI, bQ] = synthNoise(s, 24e6, round(s * 1e4) + 1, 127.5, 1);
        else
            [bI, bQ] = synthNoise(s, 24e6, round(s * 1e4) + 2, 127.397, 1.05);
        end

        H = byteHist(bI, bQ);
        sRef = s;
        if variant == "imbalance", sRef = s * sqrt((1 + 1.05^2) / 2); end
        [full, names] = allEstimators(H);

        for e = 1:numel(names)
            row = table(variant, s, string(names{e}), full(e), ...
                        20 * log10(full(e) / sRef), NaN, NaN, ...
                        'VariableNames', {'variant', 'sigma', 'estimator', 'estFull', ...
                                          'biasDb', 'blockBiasDb', 'blockStdDb'});
            E = [E; row];
        end
    end
end

% Then we save both tables and print whether the two checks passed.

writetable(E, fullfile(cfg.outDir, 'synth_est.csv'));
writetable(D, fullfile(cfg.outDir, 'synth_dump.csv'));

fprintf('\nGate synthetic cell fraction = model (within 4 SE):      %s\n', passfail(gateCell));
fprintf('Gate p40 level = quantized-model level:                 %s\n', passfail(gateP40));

% Lastly, we print each estimator's error and spread at every sigma.

grid = [0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.5 2 2.5 3 4 5 6];

fprintf('\nBias of the whole-file estimate, dB (20 log10 est/true), base variant\n');
printPivot(E(E.variant == "base", :), grid, names, 'biasDb');

fprintf('\nSpread across 1 s blocks, dB (std), base variant\n');
printPivot(E(E.variant == "base", :), grid, names, 'blockStdDb');

fprintf('\nBias, centred DC variant (mu = 127.5)\n');
printPivot(E(E.variant == "centred", :), [0.6 0.8 1.0 1.5 2 4 6], names, 'biasDb');

fprintf('\nBias, 5 %% I/Q imbalance variant\n');
printPivot(E(E.variant == "imbalance", :), [0.6 0.8 1.0 1.5 2 4 6], names, 'biasDb');

fprintf('\ndump1090 on the base twins\n');
disp(D);

function t = passfail(b)

% Here we turn true or false into PASS or FAIL.

if b, t = 'PASS'; else, t = 'FAIL'; end

end

function printPivot(E, grid, names, field)

% Here we print one row per estimator and one column per sigma.

if strcmp(field, 'blockStdDb'), fmt = '%8.3f'; else, fmt = '%8.2f'; end
fprintf('%-12s', 'sigma');
fprintf('%8.2f', grid);  fprintf('\n');

for e = 1:numel(names)
    fprintf('%-12s', names{e});

    for s = grid
        r = E(abs(E.sigma - s) < 1e-9 & E.estimator == names{e}, :);

        if isempty(r) || isnan(r.(field)(1)), fprintf('%8s', '-');
        else, fprintf(fmt, r.(field)(1)); end

    end

    fprintf('\n');

end

end
