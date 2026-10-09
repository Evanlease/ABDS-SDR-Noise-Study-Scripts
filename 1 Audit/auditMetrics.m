% First we load dump1090's strength table and the step 1 results.

cfg = expConfig();
lut = dumpLUT();  lut2 = lut.^2;
n = numel(cfg.names);
S1 = readtable(fullfile(cfg.outDir, 'step1.csv'));

% Then we create an empty table for the results.

vars = {'capture','reported','totalDb','cellPct','sigmaDC','sigmaI','sigmaQ', ...
        'bgModelDb','inflationDb','modelResidDb','p40All','p40RangeDb', ...
        'p40BlockMin','p40BlockMax','rangeFinalDb','cumInnerPct','keptPct'};

R = table('Size', [n numel(vars)], ...
          'VariableTypes', [{'string'}, repmat({'double'}, 1, numel(vars)-1)], ...
          'VariableNames', vars);

% Start a for loop over the captures, where we read each file and count the
% samples in the 4 center squares.

for c = 1:n
    name = cfg.names{c};
    f = fullfile(cfg.iqDir, [name '.iq']);
    [bI, bQ] = readUC8(f);

    H = byteHist(bI, bQ);
    N = sum(H(:));
    muI = S1.meanI(c);  muQ = S1.meanQ(c);
    P = sum(sum(H(128:129, 128:129))) / N;

    % Then we find sigma from that count, for I and Q together and for each one
    % alone.

    s = sigmaFromCellDC(P, muI, muQ);
    PI = sum(sum(H(128:129, :))) / N;  PQ = sum(sum(H(:, 128:129))) / N;
    sI = sigmaFromCellDC(PI^2, muI, muI);
    sQ = sigmaFromCellDC(PQ^2, muQ, muQ);

    % Next we find the background power that sigma should give.

    pI = quantModel(muI, s);  pQ = quantModel(muQ, s);
    bgPow = sum(sum((pI * pQ') .* lut2)) / 65535^2;

    % After that we find what dump1090's automatic gain would read.

    T = readDecoded(name);
    A = adaptiveP40(bI, bQ, T, lut);

    % Finish off the loop by filling in this capture's row.

    R.capture(c) = cfg.labels{c};
    R.reported(c) = S1.noiseReported(c);
    R.totalDb(c) = S1.totalDb(c);
    R.cellPct(c) = 100 * P;

    R.sigmaDC(c) = s;  R.sigmaI(c) = sI;  R.sigmaQ(c) = sQ;
    R.bgModelDb(c) = 10 * log10(bgPow);
    R.inflationDb(c) = S1.noiseMinusModeS(c) - R.bgModelDb(c);
    R.modelResidDb(c) = S1.totalDb(c) - R.bgModelDb(c);

    R.p40All(c) = A.levelAll;
    R.p40RangeDb(c) = -20 * log10(A.levelAll / 65536);
    R.p40BlockMin(c) = min(A.blockLevel);  R.p40BlockMax(c) = max(A.blockLevel);

    R.rangeFinalDb(c) = A.rangeDb(end);
    R.cumInnerPct(c) = 100 * A.cumInner;
    R.keptPct(c) = 100 * A.keptFraction;

    fprintf('%s done\n', cfg.labels{c});

end

writetable(R, fullfile(cfg.outDir, 'audit.csv'));

% Next we print the reported noise against the background, and I against Q.

fprintf('\n1. Reported noise vs model background (dBFS, dump1090 convention)\n');
disp(R(:, {'capture','reported','totalDb','sigmaDC','bgModelDb','inflationDb','modelResidDb'}));

fprintf('I/Q imbalance, per-component sigma from marginals:\n');
disp(R(:, {'capture','sigmaI','sigmaQ'}));

% Then we compare captures 5 and 7 with control 6, which was taken between them.

fprintf('3. Triad, control 006 as the background reference (no inversion)\n');
fprintf('   005 minus 006: %+.2f dB, 007 minus 006: %+.2f dB (reported, unrounded)\n', ...
        S1.noiseMinusModeS(5) - S1.noiseMinusModeS(6), ...
        S1.noiseMinusModeS(7) - S1.noiseMinusModeS(6));
fprintf('   cell fractions 005/006/007: %.2f / %.2f / %.2f %%\n', R.cellPct(5:7));

% After that we print what the automatic gain's 40th percentile reads on each
% capture.

fprintf('\n4. Adaptive gain p40 (magnitude levels of 65535; 290/369/434 are the three innermost)\n');
disp(R(:, {'capture','p40All','p40RangeDb','p40BlockMin','p40BlockMax','rangeFinalDb','cumInnerPct','keptPct'}));

% Then we start a for loop that steps sigma through the model to find where that
% 40th percentile first changes.

mu = mean([S1.meanI; S1.meanQ]);
sg = 0.40:0.001:1.50;
lev = zeros(size(sg));

for k = 1:numel(sg)
    p = quantModel(mu, sg(k));
    J = p * p';

    [ml, ~, g] = unique(lut(:));
    pm = accumarray(g, J(:));
    lev(k) = ml(find(cumsum(pm) > 0.40, 1));

end

chg = find(diff(lev) ~= 0);
fprintf('Model p40 transitions (mu = %.4f):\n', mu);

for k = chg
    fprintf('   sigma %.3f -> %.3f LSB: level %d -> %d (range %.2f -> %.2f dB)\n', ...
            sg(k), sg(k+1), lev(k), lev(k+1), ...
            -20*log10(lev(k)/65536), -20*log10(lev(k+1)/65536));
    if sg(k) > 1.2, break; end
end

% Lastly, we work out what background PiAware's 30 dB target aims for.

lvl30 = 65536 * 10^(-30/20);
fprintf('Target 30 dB: p40 = %.0f of 65535 = %.2f LSB magnitude -> sigma = %.2f LSB (Rayleigh p40 = 1.0108 sigma)\n', ...
        lvl30, lvl30 / 65535 * 128, lvl30 / 65535 * 128 / sqrt(-2*log(0.6)));
