% First we load earlier results and set up the thresholds, the expected counts,
% and the guard widths (extra samples to hide around each message).

cfg = expConfig();
aud = readtable(fullfile(cfg.outDir, 'audit.csv'));
acf = readtable(fullfile(cfg.outDir, 'real_acf.csv'));

ctl = ismember(string(acf.capture), ["004 control", "006 control"]);
target = mean([acf.II1(ctl), acf.II2(ctl), acf.II3(ctl)], 1);

lut = dumpLUT();
Pt = (lut / 65535).^2;

T4 = [3 5 10 20];
ledger = [327648 151344 42834 4170; 329394 167169 34199 207; 255262 94150 26075 1534; ...
          9603 0 0 0; 204513 113941 55959 29697; 7628 0 0 0; 143862 97157 64481 29268];

guards = [0 10 25 50 100 250 500];
dataLen = containers.Map({'long', 'short'}, {floor(112 * 12 / 5), floor(56 * 12 / 5)});

% Start a for loop over the captures, where we read each file and find the power
% of every sample.

Nn = table();  Tt = table();  gateOK = true;

for c = 1:numel(cfg.names)
    [bI, bQ] = readUC8(fullfile(cfg.iqDir, [cfg.names{c} '.iq']), 1, cfg.nSamples);
    N = numel(bI);

    m2 = (double(bI) - 127.5).^2 + (double(bQ) - 127.5).^2;
    P = Pt(sub2ind([256 256], double(bI) + 1, double(bQ) + 1));
    rmsAll = sqrt(mean(m2));

    % Then we check that with nothing hidden we get the expected counts of loud
    % samples.

    got = arrayfun(@(t) sum(m2 > (t * rmsAll)^2), T4);
    ok = isequal(got, ledger(c, :));  gateOK = gateOK && ok;
    fprintf('%s  empty-mask counts %s  ledger %s  %s\n', cfg.labels{c}, mat2str(got), ...
            mat2str(ledger(c, :)), passFail(ok));

    % Next we find the background level from sigma and build a twin of made-up
    % noise to compare.

    s = aud.sigmaDC(c);
    p = quantModel(127.397, s);
    rmsBg = sqrt(2 * sum(p .* ((0:255)' - 127.5).^2));
    bgDb = 10 * log10(sum(sum((p * p') .* Pt)));

    [tI, tQ] = colorTwin(s, N, round(s * 1e4) + 3, target);
    t2 = (double(tI) - 127.5).^2 + (double(tQ) - 127.5).^2;
    clear tI tQ

    % After that we copy dump1090's own noise number: total power minus decoded
    % Mode S messages.

    D = readDecoded(cfg.names{c});
    isS = D.kind ~= "modeac";
    sub = 0;

    for k = find(isS)'
        j0 = D.sample_start(k) + 1 + cfg.preambleSamples;
        idx = j0 : min(j0 + dataLen(char(D.kind(k))) - 1, N);
        sub = sub + sum(P(idx));
    end

    emuDb = 10 * log10((sum(P) - sub) / N);

    % Finish off the loop by hiding decoded messages with wider and wider
    % guards, and finding how much power and how many loud samples are left.

    for set = ["modeS", "modeS+AC"]
        use = isS;  if set == "modeS+AC", use = true(height(D), 1); end

        for g = guards
            keep = true(N, 1);

            for k = find(use)'
                keep(max(D.sample_start(k) + 1 - g, 1) : min(D.sample_end(k) + g, N)) = false;
            end

            nk = sum(keep);
            keptDb = 10 * log10(mean(P(keep)));

            Nn = [Nn; table(string(cfg.labels{c}), set, g, 1 - nk / N, keptDb, bgDb, keptDb - bgDb, ...
                            emuDb, aud.reported(c), ...
                            'VariableNames', {'capture', 'maskSet', 'guard', 'maskedFrac', 'keptDbfs', ...
                                              'bgModelDbfs', 'leakDb', 'dumpEmuDbfs', 'reportedDbfs'})];

            for conv = ["ledger", "background"]
                if conv == "ledger", sc = rmsAll; else, sc = rmsBg; end

                for t = T4
                    thr = (t * sc)^2;
                    n = sum(m2(keep) > thr);
                    nul = sum(t2 > thr) * nk / N;

                    Tt = [Tt; table(string(cfg.labels{c}), set, g, conv, t, sc, n, nul, n / max(nul, 1), ...
                                    'VariableNames', {'capture', 'maskSet', 'guard', 'convention', 't', ...
                                                      'scale', 'count', 'twinCount', 'ratio'})];

                end
            end
        end
    end

    fprintf('%s  rmsAll %.4f  rmsBg %.4f  dump emu %.2f vs reported %.1f dBFS\n', ...
            cfg.labels{c}, rmsAll, rmsBg, emuDb, aud.reported(c));

end

% After the loop we print the check result and save both tables.

fprintf('Gate, empty mask reproduces the ledger: %s\n', passFail(gateOK));

writetable(Nn, fullfile(cfg.outDir, 'mask_noise.csv'));
writetable(Tt, fullfile(cfg.outDir, 'mask_tail.csv'));

% Lastly, we print the leftover power for each guard width, and the loud samples
% at 3x.

labs = string(cfg.labels);
fprintf('\nLeak = kept mean power - background model, dB (dump1090 convention)\n');

for set = ["modeS", "modeS+AC"]
    fprintf('%s\n%-14s', set, 'guard'); fprintf('%8d', guards); fprintf('\n');

    for c = 1:numel(labs)
        fprintf('%-14s', labs(c));
        fprintf('%8.3f', Nn.leakDb(Nn.capture == labs(c) & Nn.maskSet == set));
        fprintf('\n');
    end
end

fprintf('\nTail at 3x, modeS+AC mask: count / twin null, by guard\n');

for conv = ["ledger", "background"]
    fprintf('%s\n%-14s', conv, 'guard'); fprintf('%8d', guards); fprintf('\n');

    for c = 1:numel(labs)
        fprintf('%-14s', labs(c));
        fprintf('%8.2f', Tt.ratio(Tt.capture == labs(c) & Tt.maskSet == "modeS+AC" & ...
                                  Tt.convention == conv & Tt.t == 3));
        fprintf('\n');
    end
end

function s = passFail(ok)

% Here we turn true or false into PASS or FAIL.

if ok, s = 'PASS'; else, s = 'FAIL'; end

end
