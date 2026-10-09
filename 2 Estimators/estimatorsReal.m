% First we set up 30 blocks of 1 s each and two empty tables.

cfg = expConfig();
nBlk = 30;  blkLen = 2.4e6;
R = table();  A = table();

% Start a for loop over the captures, where we run every estimator on the whole
% file.

for c = 1:numel(cfg.names)
    [bI, bQ] = readUC8(fullfile(cfg.iqDir, [cfg.names{c} '.iq']), 1, cfg.nSamples);
    N = numel(bI);
    [full, names] = allEstimators(byteHist(bI, bQ));

    % Then we run them again with every decoded message cut out.

    T = readDecoded(cfg.names{c});
    keep = true(N, 1);

    for k = 1:height(T)
        keep(T.sample_start(k) + 1 : min(T.sample_end(k), N)) = false;
    end

    masked = allEstimators(byteHist(bI(keep), bQ(keep)));

    % Next we run them on each 1 s block.

    blk = nan(nBlk, numel(names));

    for k = 1:nBlk
        idx = (k - 1) * blkLen + (1:blkLen);
        blk(k, :) = allEstimators(byteHist(bI(idx), bQ(idx)));
    end

    % Finish off the loop by adding a row per estimator, plus how alike nearby
    % samples are (the autocorrelation).

    for e = 1:numel(names)
        R = [R; table(string(cfg.labels{c}), cfg.is1090(c), string(names{e}), ...
                      full(e), masked(e), 20 * log10(full(e) / masked(e)), ...
                      std(20 * log10(blk(:, e))), ...
                      20 * log10(max(blk(:, e)) / min(blk(:, e))), ...
                      'VariableNames', {'capture', 'is1090', 'estimator', 'estFull', ...
                                        'estMasked', 'fullMinusMaskedDb', ...
                                        'blockStdDb', 'blockRangeDb'})];
    end

    A = [A; acfRow(string(cfg.labels{c}), bI, bQ)];
    fprintf('%s  masked out %.2f %%\n', cfg.labels{c}, 100 * mean(~keep));

end

% Next we add the autocorrelation of a made-up noise twin of capture 6.

aud = readtable(fullfile(cfg.outDir, 'audit.csv'));
s6 = round(aud.sigmaDC(6), 4);
[bI, bQ] = synthNoise(s6, cfg.nSamples, round(s6 * 1e4));
A = [A; acfRow("twin of 006", bI, bQ)];

% Then we save both tables.

writetable(R, fullfile(cfg.outDir, 'real_est.csv'));
writetable(A, fullfile(cfg.outDir, 'real_acf.csv'));

% After that we print every estimator's sigma for each capture.

labs = string(cfg.labels);
fprintf('\nWhole-buffer estimate, sigma in LSB\n%-12s', '');
fprintf('%13s', labs); fprintf('\n');

for e = 1:numel(names)
    fprintf('%-12s', names{e});

    for c = 1:numel(labs)
        fprintf('%13.4f', R.estFull(R.capture == labs(c) & R.estimator == names{e}));
    end

    fprintf('\n');

end

% Then we print how much traffic moves each estimate: captures 5 and 7 minus
% control 6.

fprintf('\nTraffic bias vs the 006 control, dB (0 = immune)\n');
fprintf('%-12s%10s%10s%12s%12s\n', '', '005-006', '007-006', 'max|fm| 1090', 'blockStd 5/7');

for e = 1:numel(names)
    g = @(lab, col) R.(col)(R.capture == lab & R.estimator == names{e});
    d5 = 20 * log10(g(labs(5), 'estFull') / g(labs(6), 'estFull'));
    d7 = 20 * log10(g(labs(7), 'estFull') / g(labs(6), 'estFull'));
    fm = R.fullMinusMaskedDb(R.estimator == names{e} & R.is1090 == 1);

    fprintf('%-12s%10.3f%10.3f%12.3f%6.3f%6.3f\n', names{e}, d5, d7, max(abs(fm)), ...
            g(labs(5), 'blockStdDb'), g(labs(7), 'blockStdDb'));

end

% Lastly, we print how much cutting out decoded messages changes each estimate.

fprintf('\nFull minus masked (decoded traffic effect), dB\n%-12s', '');
fprintf('%13s', labs); fprintf('\n');

for e = 1:numel(names)
    fprintf('%-12s', names{e});

    for c = 1:numel(labs)
        fprintf('%13.3f', R.fullMinusMaskedDb(R.capture == labs(c) & R.estimator == names{e}));
    end

    fprintf('\n');

end

% We also print the autocorrelation table.

fprintf('\nByte autocorrelation\n');
disp(A)

function row = acfRow(lab, bI, bQ)

% Here we find the autocorrelation 1 to 3 samples apart, plus I against Q.

x = double(bI) - mean(double(bI));  y = double(bQ) - mean(double(bQ));
vx = mean(x.^2);  vy = mean(y.^2);

r = @(a, b, L, va, vb) mean(a(1:end-L) .* b(1+L:end)) / sqrt(va * vb);
row = table(lab, r(x, x, 1, vx, vx), r(x, x, 2, vx, vx), r(x, x, 3, vx, vx), ...
            r(y, y, 1, vy, vy), r(x, y, 0, vx, vy), r(x, y, 1, vx, vy), ...
            'VariableNames', {'capture', 'II1', 'II2', 'II3', 'QQ1', 'IQ0', 'IQ1'});

end
