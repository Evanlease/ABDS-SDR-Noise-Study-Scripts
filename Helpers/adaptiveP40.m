function A = adaptiveP40(bI, bQ, T, lut)

% First we set up dump1090's automatic gain settings.

cfg = expConfig();
N = numel(bI);
pct = 40;  alpha = 1/3;
blockLen = round(cfg.fs);

% Then we create a mask that hides the samples of every decoded Mode S message,
% like dump1090 does.

keep = true(N, 1);

for k = 1:height(T)
    if T.kind(k) == "modeac", continue; end
    if T.kind(k) == "long", bits = 112; else, bits = 56; end

    a = T.sample_start(k) + 1;
    b = min(T.sample_start(k) + floor((bits + 8) * 12 / 5), N);
    keep(a:b) = false;

end

% Next we turn each byte pair into dump1090's signal strength.

mag = lut(sub2ind([256 256], double(bI) + 1, double(bQ) + 1));

% Then we start a for loop that finds the 40th percentile in each 1 s block.

nBlocks = floor(N / blockLen);
A.blockLevel = zeros(nBlocks, 1);

for blk = 1:nBlocks
    idx = (blk - 1) * blockLen + (1:blockLen);
    m = mag(idx(keep(idx)));
    A.blockLevel(blk) = pctLevel(m, pct);
end

% After that we smooth the block levels and turn them into the range in dB that
% the gain control uses.

A.smoothed = zeros(nBlocks, 1);
A.smoothed(1) = A.blockLevel(1);

for blk = 2:nBlocks
    A.smoothed(blk) = A.smoothed(blk-1) * (1 - alpha) + A.blockLevel(blk) * alpha;
end

A.rangeDb = -20 * log10(A.smoothed / 65536);

% Lastly, we get the 40th percentile of the whole capture and how much of it we
% kept.

mk = mag(keep);
A.levelAll = pctLevel(mk, pct);
A.cumInner = mean(mk == lut(128, 128) | mk == lut(128, 129) | mk == lut(129, 128));
A.keptFraction = mean(keep);

end

function level = pctLevel(m, pct)

% Here we find the level that pct percent of the samples sit below.

counts = accumarray(m + 1, 1, [65536 1]);
countN = floor(numel(m) * pct / 100);
level = find(cumsum(counts) > countN, 1) - 1;

end
