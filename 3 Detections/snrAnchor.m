% First we load earlier results and build the receiver filter.

cfg = expConfig();
acf = readtable(fullfile(cfg.outDir, 'real_acf.csv'));
rest = readtable(fullfile(cfg.outDir, 'real_est.csv'));
inj = readtable(fullfile(cfg.outDir, 'injection.csv'));

ctl = ismember(string(acf.capture), ["004 control", "006 control"]);
target = mean([acf.II1(ctl), acf.II2(ctl), acf.II3(ctl)], 1);
sigCF = rest.estFull(string(rest.estimator) == "centralFit");
[~, ~, h] = colorTwin(sigCF(6), 10, 1, target);

mu = 127.397;  lead = 480;

% Then we calibrate: we add made-up messages of known SNR (signal-to-noise
% ratio) to noise, let dump1090 find them, and measure each one's SNR ourselves.

s = 0.72;  N = 12e6;  lv = [4 6 8 10 12 14 17 20 25 30];  nPer = 300;
rng(2468, 'twister');

lev = repmat(1:numel(lv), 1, nPer);  lev = lev(randperm(numel(lev)))';
F = numel(lev);  starts = 5000 + (0:F-1)' * 3400;

frac = rand(F, 1);  phi = 2 * pi * rand(F, 1);
bits = false(F, 112);  hexes = strings(F, 1);

for f = 1:F
    [bits(f, :), hexes(f)] = modesFrame(hex2dec('AC0000') + lev(f), sprintf('%02d%06d', lev(f), f));
end

gI = s * filter(h, 1, randn(N + 3, 1));  gI = gI(4:end);
gQ = s * filter(h, 1, randn(N + 3, 1));  gQ = gQ(4:end);

[gI, gQ] = injectFrames(gI, gQ, starts, bits, s * sqrt(2 * 10 .^ (lv(lev)' / 10)), phi, frac, h);
bI = uint8(min(max(round(mu + gI), 0), 255));  bQ = uint8(min(max(round(mu + gQ), 0), 255));
clear gI gQ

writeUC8(cfg.work, bI, bQ);
msg = runDump1090(cfg.work, 'raw', '--fix');
[tf, loc] = ismember(msg.hex, hexes);

zI = double(bI) - mu;  zQ = double(bQ) - mu;
C = table();

for k = find(tf)'
    f = loc(k);
    st = floor(msg.ts(k) / 5) - lead;
    [snr, tau] = pulseSnr(zI, zQ, st, bits(f, :), s, h);

    C = [C; table(lv(lev(f)), snr, st + tau - (starts(f) + frac(f)), ...
                  'VariableNames', {'trueSnrDb', 'estSnrDb', 'timingErr'})];

end

C = unique(C, 'rows', 'stable');
writetable(C, fullfile(cfg.outDir, 'snr_anchor_cal.csv'));

% Next we print how far our measured SNR is from the true SNR at each level.

fprintf('Calibration, sigma %.2f, decoded synthetic frames\n', s);
fprintf('%8s %6s %10s %10s %10s %12s\n', 'true dB', 'n', 'median est', 'bias', 'IQR', 'median |dt|');

for v = lv
    r = C(C.trueSnrDb == v, :);
    if isempty(r), continue; end
    fprintf('%8d %6d %10.2f %10.2f %10.2f %12.2f\n', v, height(r), median(r.estSnrDb), ...
            median(r.estSnrDb) - v, iqr(r.estSnrDb), median(abs(r.timingErr)));
end

% Then we start a for loop that measures the SNR of every real Mode S message
% dump1090 decoded.

R = table();

for c = find(cfg.is1090)
    [bI, bQ] = readUC8(fullfile(cfg.iqDir, [cfg.names{c} '.iq']), 1, cfg.nSamples);
    zI = double(bI) - mu;  zQ = double(bQ) - mu;

    D = readDecoded(cfg.names{c});
    D = D(D.kind ~= "modeac", :);

    for k = 1:height(D)
        b = hex2bits(D.hex(k));
        snr = pulseSnr(zI, zQ, D.sample_start(k), b, sigCF(c), h);
        R = [R; table(string(cfg.labels{c}), D.kind(k), D.df(k), snr, ...
                      'VariableNames', {'capture', 'kind', 'df', 'snrDb'})];
    end
end

writetable(R, fullfile(cfg.outDir, 'snr_anchor_real.csv'));

% Lastly, we print the SNR of the weakest real messages and compare with the
% injection.

pc = [1 5 10 25 50];
fprintf('\nReal decoded Mode S, estimated SNR (dB) percentiles %s\n', mat2str(pc));

for lab = unique(R.capture, 'stable')'
    for kd = ["long", "short"]
        x = R.snrDb(R.capture == lab & R.kind == kd);
        fprintf('%-14s %-5s n=%5d  %s\n', lab, kd, numel(x), sprintf('%7.2f', prctile(x, pc)));
    end
end

for kd = ["long", "short"]
    x = R.snrDb(ismember(R.capture, ["001 morning", "002 midday", "003 afternoon"]) & R.kind == kd);
    fprintf('001-003 pooled %-5s n=%5d  %s\n', kd, numel(x), sprintf('%7.2f', prctile(x, pc)));
end

q = inj(inj.source == "synthetic" & abs(inj.sigma - 0.72) < 1e-9, :);
fprintf('\nInjection, sigma 0.72: P(decode) at 8/10/12/14 dB = %.3f %.3f %.3f %.3f\n', ...
        q.pDecode(q.snrDb == 8), q.pDecode(q.snrDb == 10), q.pDecode(q.snrDb == 12), q.pDecode(q.snrDb == 14));

function [snrDb, tau] = pulseSnr(zI, zQ, st, bits, sigma, h)

% Here we cut out the power around the message and set up when each pulse
% starts.

nb = numel(bits);
t0 = [0 1.0 3.5 4.5, 8 + (0:nb-1) + 0.5 * ~bits];
L = ceil((8 + nb) * 2.4) + 2 + numel(h);
pad = 4;

idx = st - pad + (0:L + 2*pad - 1)' + 1;
idx = idx(idx >= 1 & idx <= numel(zI));
p = zI(idx).^2 + zQ(idx).^2;
n = (idx - 1 - st);

% Then we start a for loop that slides the expected pulse pattern over the
% samples, keeps the best fit, and turns it into SNR.

best = -Inf;  tau = 0;  A2 = 0;
grid = -3:0.25:3;

for pass = 1:2
    for tt = grid
        w = template(t0, tt, n, h);
        off = w.^2 < 1e-6;
        if sum(off) < 10, continue; end

        pOff = mean(p(off));
        num = sum(w.^2 .* (p - pOff));  den = sum(w.^4);
        score = num^2 / den * sign(num);
        if score > best, best = score; tau = tt; A2 = num / den; end

    end

    grid = tau + (-0.25:0.05:0.25);

end

snrDb = 10 * log10(max(A2, 1e-6) / (2 * sigma^2));

end

function w = template(t0, tau, n, h)

% Here we find how much of each sample the pulses cover, then filter it.

a = t0 * 2.4 + tau;  b = a + 1.2;
w = sum(max(0, min(n + 1, b) - max(n, a)), 2);
w = filter(h, 1, w);

end

function b = hex2bits(hx)

% Here we turn each hex digit into its 4 bits.

hx = char(hx);
b = false(1, 4 * numel(hx));

for k = 1:numel(hx)
    b(4*k-3:4*k) = bitget(hex2dec(hx(k)), 4:-1:1) == 1;
end

end
