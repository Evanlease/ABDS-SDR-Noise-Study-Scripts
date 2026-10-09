% First we load earlier results and build the receiver filter.

cfg = expConfig();
acf = readtable(fullfile(cfg.outDir, 'real_acf.csv'));
rest = readtable(fullfile(cfg.outDir, 'real_est.csv'));

ctl = ismember(string(acf.capture), ["004 control", "006 control"]);
target = mean([acf.II1(ctl), acf.II2(ctl), acf.II3(ctl)], 1);
sigCF = rest.estFull(string(rest.estimator) == "centralFit");

[~, ~, h] = colorTwin(sigCF(6), 10, 1, target);
lut = dumpLUT();

% Then we set up the test: 10 s files, three noise levels, and messages filling
% 0 to 20%.

N = 24e6;
mu = 127.397;  L = ceil(120 * 2.4) + 2 + numel(h);

sigmas = [0.72 1.5 4];
occ = [0 0.5 1 2 5 10 20];
pick = [1 2 3 6 8];

% Next we create the test messages and their start times, which every noise
% level shares.

rng(4242, 'twister');
Fmax = round(max(occ) / 100 * N / 288);
pool = false(Fmax, 112);

for f = 1:Fmax
    cs = char(randi([65 90], 1, 8));
    pool(f, :) = modesFrame(randi([1, 2^24 - 1]), cs);
end

sched = cell(numel(occ), 1);

for o = 1:numel(occ)
    F = round(occ(o) / 100 * N / 288);

    sched{o}.starts = sort(randi([0, N - L - 1], F, 1));
    sched{o}.snr = 3 + 19 * rand(F, 1);
    sched{o}.phi = 2 * pi * rand(F, 1);
    sched{o}.frac = rand(F, 1);

end

% Then we start a for loop where, for each noise level and fill, we add the
% messages and run the estimators and dump1090.

R = table();

for s = sigmas
    rng(round(s * 1e4) + 7, 'twister');
    g0I = s * filter(h, 1, randn(N + 3, 1));  g0I = g0I(4:end);
    g0Q = s * filter(h, 1, randn(N + 3, 1));  g0Q = g0Q(4:end);

    pI = quantModel(mu, s);  pQ = quantModel(mu, s);
    modelDb = 10 * log10(sum(sum((pI * pQ') .* lut.^2)) / 65535^2);

    for o = 1:numel(occ)
        S = sched{o};  F = numel(S.starts);
        A = s * sqrt(2 * 10 .^ (S.snr / 10));

        [gI, gQ] = injectFrames(g0I, g0Q, S.starts, pool(1:F, :), A, S.phi, S.frac, h);
        bI = uint8(min(max(round(mu + gI), 0), 255));
        bQ = uint8(min(max(round(mu + gQ), 0), 255));
        clear gI gQ

        covered = false(N, 1);

        for f = 1:F
            covered(S.starts(f) + (1:L)) = true;
        end

        est = allEstimators(byteHist(bI, bQ));
        estOut = allEstimators(byteHist(bI(~covered), bQ(~covered)));

        writeUC8(cfg.work, bI, bQ);
        st = runDump1090(cfg.work, 'stats', '--fix');

        b = 20 * log10(est(pick) / s);
        R = [R; table(s, occ(o), F, mean(covered), b(1), b(2), b(3), b(4), b(5), ...
                      20 * log10(estOut(8) / s), st.noiseDbfs, modelDb, st.accepted, st.repaired, st.preambles, ...
                      'VariableNames', {'sigma', 'occPct', 'frames', 'coveredFrac', ...
                                        'meanPowerDb', 'sheppardDb', 'p40Db', 'cellDCDb', 'centralFitDb', ...
                                        'centralFitOracleDb', 'dumpNoiseDbfs', 'modelNoiseDbfs', ...
                                        'accepted', 'repaired', 'preambles'})];

        fprintf('sigma %.2f occ %4.1f %%  cellDC %+.3f  centralFit %+.3f  meanPower %+.3f  dump %.2f dBFS\n', ...
                s, occ(o), b(4), b(5), b(1), st.noiseDbfs);

    end
end

% After that we find how much dump1090's noise rose and the share of messages it
% decoded.

R.dumpExcessDb = nan(height(R), 1);

for s = sigmas
    k = R.sigma == s;
    R.dumpExcessDb(k) = R.dumpNoiseDbfs(k) - R.dumpNoiseDbfs(k & R.occPct == 0);
end

R.acceptedFrac = (R.accepted + R.repaired) ./ R.frames;
R.acceptedFrac(R.frames == 0) = NaN;

% Lastly, we save and print the table.

writetable(R, fullfile(cfg.outDir, 'occupancy.csv'));

fprintf('\nBias vs true sigma, dB (as sigma; power bias is the same number, 20 log10 on sigma)\n');
disp(R(:, {'sigma', 'occPct', 'coveredFrac', 'cellDCDb', 'centralFitDb', 'meanPowerDb', ...
           'sheppardDb', 'p40Db', 'dumpExcessDb', 'acceptedFrac'}))
