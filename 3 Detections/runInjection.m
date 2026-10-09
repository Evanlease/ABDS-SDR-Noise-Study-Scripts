% First we load earlier results and build the receiver filter.

cfg = expConfig();
acf = readtable(fullfile(cfg.outDir, 'real_acf.csv'));
rest = readtable(fullfile(cfg.outDir, 'real_est.csv'));

ctl = ismember(string(acf.capture), ["004 control", "006 control"]);
target = mean([acf.II1(ctl), acf.II2(ctl), acf.II3(ctl)], 1);
sigCF = rest.estFull(string(rest.estimator) == "centralFit");
[~, ~, h] = colorTwin(sigCF(6), 10, 1, target);

N = cfg.nSamples;  mu = 127.397;  icaoBase = hex2dec('AB0000');

% Then we plan the test messages: 23 SNR (signal-to-noise ratio) levels, 900
% messages each.

snrDb = 0:22;  nL = numel(snrDb);  nPer = 900;  spacing = 3400;
rng(12345, 'twister');

lev = repmat(1:nL, 1, nPer);  lev = lev(randperm(numel(lev)))';
F = numel(lev);

starts = 5000 + (0:F-1)' * spacing;
frac = rand(F, 1);  phi = 2 * pi * rand(F, 1);
bits = false(F, 112);  hexes = strings(F, 1);

for f = 1:F
    [bits(f, :), hexes(f)] = modesFrame(icaoBase + lev(f), sprintf('%02d%06d', lev(f), f));
end

% Next we run the checks: strong messages must decode, and no test message can
% show up in the real captures.

[bI, bQ] = gateFile(1.0, 30, 2000, 99, [], icaoBase, h, mu, N);
writeUC8(cfg.work, bI, bQ);
[pGate, ~] = matchRate(runDump1090(cfg.work, 'raw', '--fix'), cfg, 99, icaoBase);

[bI, bQ, flipHex] = gateFile(1.0, 30, 1000, 98, 'flip', icaoBase, h, mu, N);
writeUC8(cfg.work, bI, bQ);
m0 = runDump1090(cfg.work, 'raw');  m1 = runDump1090(cfg.work, 'raw', '--fix');
fixDefault = mean(ismember(flipHex, m0.hex));  fixOn = mean(ismember(flipHex, m1.hex));

clean = 0;

for c = 1:numel(cfg.names)
    m = runDump1090(fullfile(cfg.iqDir, [cfg.names{c} '.iq']), 'raw', '--fix');
    clean = clean + sum(startsWith(m.hex, "8DAB00") | startsWith(m.hex, "8CAB00"));
end

fprintf('Gate 30 dB decode at sigma 1: %.4f (need >= 0.999): %s\n', pGate, passFail(pGate >= 0.999));
fprintf('Gate zero test-ICAO in the seven real captures: %d: %s\n', clean, passFail(clean == 0));
fprintf('Single-bit-error frames decoded: default flags %.3f, with --fix %.3f\n', fixDefault, fixOn);

% Then we start a for loop where, for each noise level, we add the messages and
% record the share dump1090 decodes.

sigmas = [0.5 0.6 0.72 0.85 1.0 1.2 1.5 2 3 4 6];
R = table();

for s = sigmas
    rng(round(s * 1e4) + 99, 'twister');
    gI = s * filter(h, 1, randn(N + 3, 1));  gI = gI(4:end);
    gQ = s * filter(h, 1, randn(N + 3, 1));  gQ = gQ(4:end);

    A = s * sqrt(2 * 10 .^ (snrDb(lev)' / 10));
    [gI, gQ] = injectFrames(gI, gQ, starts, bits, A, phi, frac, h);
    clipped = mean(abs(gI + mu - 127.5) > 127.5 | abs(gQ + mu - 127.5) > 127.5);

    bI = uint8(min(max(round(mu + gI), 0), 255));  bQ = uint8(min(max(round(mu + gQ), 0), 255));
    clear gI gQ

    writeUC8(cfg.work, bI, bQ);
    p = levelRates(runDump1090(cfg.work, 'raw', '--fix'), hexes, lev, nL, nPer);
    R = [R; table(repmat("synthetic", nL, 1), repmat(s, nL, 1), snrDb', p, ...
                  repmat(clipped, nL, 1), ...
                  'VariableNames', {'source', 'sigma', 'snrDb', 'pDecode', 'clippedFrac'})];

    fprintf('sigma %.2f  P(decode) at 6/10/14 dB: %.3f %.3f %.3f\n', s, p(7), p(11), p(15));

end

% After that we put the same messages into real capture 6 and its twin, as a
% check.

s6 = sigCF(6);
[rI, rQ] = readUC8(fullfile(cfg.iqDir, [cfg.names{6} '.iq']), 1, N);
[tI, tQ] = colorTwin(s6, N, round(s6 * 1e4), target);
A = s6 * sqrt(2 * 10 .^ (snrDb(lev)' / 10));

for src = ["real 006", "twin 006"]
    if src == "real 006", gI = double(rI) - mu; gQ = double(rQ) - mu;
    else, gI = double(tI) - mu; gQ = double(tQ) - mu; end

    [gI, gQ] = injectFrames(gI, gQ, starts, bits, A, phi, frac, h);
    bI = uint8(min(max(round(mu + gI), 0), 255));  bQ = uint8(min(max(round(mu + gQ), 0), 255));

    writeUC8(cfg.work, bI, bQ);
    p = levelRates(runDump1090(cfg.work, 'raw', '--fix'), hexes, lev, nL, nPer);
    R = [R; table(repmat(src, nL, 1), repmat(s6, nL, 1), snrDb', p, zeros(nL, 1), ...
                  'VariableNames', {'source', 'sigma', 'snrDb', 'pDecode', 'clippedFrac'})];

end

writetable(R, fullfile(cfg.outDir, 'injection.csv'));

% Lastly, we find the SNR where decoding reaches 50% and 90% at each noise
% level.

S = table();

for src = unique(R.source, 'stable')'
    for s = unique(R.sigma(R.source == src))'
        r = R(R.source == src & R.sigma == s, :);
        S = [S; table(src, s, crossing(r.snrDb, r.pDecode, 0.5), crossing(r.snrDb, r.pDecode, 0.9), ...
                      'VariableNames', {'source', 'sigma', 'snr50', 'snr90'})];
    end
end

ref = S.snr50(S.source == "synthetic" & S.sigma == 6);
ref90 = S.snr90(S.source == "synthetic" & S.sigma == 6);
S.loss50 = S.snr50 - ref;  S.loss90 = S.snr90 - ref90;

writetable(S, fullfile(cfg.outDir, 'injection_summary.csv'));
disp(S)

function [bI, bQ, sentHex] = gateFile(s, snr, n, levelCode, mode, icaoBase, h, mu, N)

% Here we make a test file of n messages at one SNR, with one bit flipped if
% asked.

rng(777 + levelCode, 'twister');
gI = s * filter(h, 1, randn(N + 3, 1));  gI = gI(4:end);
gQ = s * filter(h, 1, randn(N + 3, 1));  gQ = gQ(4:end);

st = 5000 + (0:n-1)' * 3400;
b = false(n, 112);  sentHex = strings(n, 1);

for f = 1:n
    [b(f, :), sentHex(f)] = modesFrame(icaoBase + levelCode, sprintf('%02d%06d', levelCode, f));

    if strcmp(mode, 'flip')
        k = 8 + randi(80);
        b(f, k) = ~b(f, k);
    end
end

A = repmat(s * sqrt(2 * 10^(snr / 10)), n, 1);
[gI, gQ] = injectFrames(gI, gQ, st, b, A, 2 * pi * rand(n, 1), rand(n, 1), h);
bI = uint8(min(max(round(mu + gI), 0), 255));  bQ = uint8(min(max(round(mu + gQ), 0), 255));

end

function [p, got] = matchRate(msg, ~, levelCode, icaoBase)

% Here we count the share of the 2000 test messages that came back.

pre = "8D" + string(dec2hex(icaoBase + levelCode, 6));
got = unique(msg.hex(startsWith(msg.hex, pre)));
p = numel(got) / 2000;

end

function p = levelRates(msg, hexes, lev, nL, nPer)

% Here we match decodes to the messages we sent and find the share decoded at
% each level.

[tf, loc] = ismember(msg.hex, hexes);
got = unique(loc(tf));
p = accumarray(lev(got), 1, [nL 1]) / nPer;

end

function x = crossing(snr, p, level)

% Here we find the SNR where the decode curve first reaches the level.

k = find(p >= level, 1);

if isempty(k), x = NaN; elseif k == 1, x = snr(1);
else, x = snr(k-1) + (level - p(k-1)) / (p(k) - p(k-1)) * (snr(k) - snr(k-1));
end

end

function s = passFail(ok)

% Here we turn true or false into PASS or FAIL.

if ok, s = 'PASS'; else, s = 'FAIL'; end

end
