% First we get dump1090's strength table and square it.

cfg = expConfig();
lut = dumpLUT();
lut2 = lut.^2;

% Then we write down the numbers we expect from earlier notes.

expPre   = [461185 400498 342520 498369 357695 359962 355392];
expBad   = [1831322 1574685 1332417 1999909 1395697 1410654 1389227];
expAcc   = [1844 1933 1591 0 1412 0 912];
expAC    = [1673 2574 1201 306 1330 12 866];

expNoise = [-39.6 -40.5 -41.2 -39.7 -40.4 -41.4 -40.5];
expRms   = [1.4834 1.3189 1.1949 1.3385 1.7983 1.1036 1.5842];
expCell  = [57.45 64.64 70.04 54.20 68.96 69.28 69.26];

% Next we set up the results table and four checks that start out passing.

n = numel(cfg.names);

R = table('Size', [n 16], ...
    'VariableTypes', [{'string'}, repmat({'double'}, 1, 15)], ...
    'VariableNames', {'capture', 'preambles', 'bad', 'accepted', 'modeac', ...
                      'noiseReported', 'meanI', 'meanQ', 'rms127p5', 'cellPct', ...
                      'totalDb', 'noiseMinusModeS', 'noiseMinusAll', ...
                      'modeSSpans', 'modeACSpans', 'spanSamples'});

ok = true(1, 4);

% Start a for loop that goes through the captures one at a time.

for c = 1:n
    name = cfg.names{c};
    f = fullfile(cfg.iqDir, [name '.iq']);
    fprintf('%s\n', cfg.labels{c});

    % For check A, we make sure dump1090's counts and noise match the notes.

    s = runDump1090(f, 'stats');

    R.capture(c) = cfg.labels{c};
    R.preambles(c) = s.preambles;   R.bad(c) = s.badFormat;
    R.accepted(c) = s.accepted;     R.modeac(c) = s.modeac;
    R.noiseReported(c) = s.noiseDbfs;

    if s.preambles ~= expPre(c) || s.badFormat ~= expBad(c) || ...
       s.accepted ~= expAcc(c) || s.modeac ~= expAC(c) || ...
       abs(s.noiseDbfs - expNoise(c)) > 0.01
        ok(1) = false;
        fprintf('  A FAIL: counters differ from DATA.md\n');
    end

    % For check B, we make sure the bytes' rms (typical size) and center share
    % match.

    [bI, bQ] = readUC8(f);
    H = byteHist(bI, bQ);
    N = sum(H(:));
    v = (0:255)';

    R.meanI(c) = sum(sum(H, 2) .* v) / N;
    R.meanQ(c) = sum(sum(H, 1)' .* v) / N;

    d = v - 127.5;
    R.rms127p5(c) = sqrt((sum(sum(H, 2) .* d.^2) + sum(sum(H, 1)' .* d.^2)) / N);
    R.cellPct(c) = 100 * sum(sum(H(128:129, 128:129))) / N;

    if abs(R.rms127p5(c) - expRms(c)) > 0.00005 || abs(R.cellPct(c) - expCell(c)) > 0.005
        ok(2) = false;
        fprintf('  B FAIL: rms %.4f vs %.4f, cell %.2f vs %.2f\n', ...
                R.rms127p5(c), expRms(c), R.cellPct(c), expCell(c));
    end

    % For check C, we rebuild dump1090's noise number as total power minus
    % decoded messages.

    totalPow = sum(H(:) .* lut2(:));
    R.totalDb(c) = 10 * log10(totalPow / N / 65535^2);

    T = readDecoded(name);
    sigS = 0; sigAC = 0; nS = 0; nAC = 0; spanN = 0;

    for k = 1:height(T)
        if T.kind(k) == "modeac"
            a = T.sample_start(k) + 1;
            b = min(T.sample_end(k), N);

            sigAC = sigAC + sum(lut2(sub2ind([256 256], ...
                        double(bI(a:b)) + 1, double(bQ(a:b)) + 1)));
            nAC = nAC + 1;

        else
            if T.kind(k) == "long", bits = 112; else, bits = 56; end
            a = T.sample_start(k) + 1 + cfg.preambleSamples;
            b = min(a + floor(bits * 12 / 5) - 1, N);

            sigS = sigS + sum(lut2(sub2ind([256 256], ...
                        double(bI(a:b)) + 1, double(bQ(a:b)) + 1)));
            nS = nS + 1;
            spanN = spanN + (b - a + 1);

        end
    end

    R.noiseMinusModeS(c) = 10 * log10((totalPow - sigS) / N / 65535^2);
    R.noiseMinusAll(c)   = 10 * log10((totalPow - sigS - sigAC) / N / 65535^2);
    R.modeSSpans(c) = nS;  R.modeACSpans(c) = nAC;  R.spanSamples(c) = spanN;

    if abs(R.noiseMinusModeS(c) - R.noiseReported(c)) > 0.05 + 1e-9
        ok(3) = false;
        fprintf('  C FAIL: rebuilt %.3f vs reported %.1f\n', ...
                R.noiseMinusModeS(c), R.noiseReported(c));
    end

    % For check D, we write capture 6 back out and make sure the copy matches
    % the original.

    if c == 6
        writeUC8(cfg.work, bI, bQ);
        st = system(sprintf('cmp -s "%s" "%s"', f, cfg.work));

        if st ~= 0
            ok(4) = false;
            fprintf('  D FAIL: round trip differs\n');
        end

        [pI, pQ] = readUC8(f, 1000001, 5000);

        if ~isequal(pI, bI(1000001:1005000)) || ~isequal(pQ, bQ(1000001:1005000))
            ok(4) = false;
            fprintf('  D FAIL: partial read differs\n');
        end
    end
end

% Lastly, we save the table and print PASS or FAIL for each check.

writetable(R, fullfile(cfg.outDir, 'step1.csv'));
disp(R(:, {'capture', 'noiseReported', 'totalDb', 'noiseMinusModeS', ...
           'noiseMinusAll', 'meanI', 'meanQ'}));

gates = {'A counters', 'B census', 'C noise rebuild', 'D round trip'};

for g = 1:4
    if ok(g), st = 'PASS'; else, st = 'FAIL'; end
    fprintf('%-16s %s\n', gates{g}, st);
end
