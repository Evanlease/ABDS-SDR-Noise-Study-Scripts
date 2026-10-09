% First we load the result tables and make sure they list the captures in the
% same order.

cfg = expConfig();

A = readtable(fullfile(cfg.outDir, 'audit.csv'),       'TextType', 'string');
F = readtable(fullfile(cfg.outDir, 'fix_compare.csv'), 'TextType', 'string');
C = readtable(fullfile(cfg.outDir, 'color_twins.csv'), 'TextType', 'string');
D = readtable(fullfile(cfg.outDir, 'dither.csv'),      'TextType', 'string');

if ~isequal(A.capture, F.capture) || ~isequal(A.capture, C.capture)
    error('audit, fix_compare and color_twins list the captures in different orders.');
end

% These are the sigmas, the real preamble counts, and the twin counts.

sigma  = A.sigmaDC;
realN  = F.preambles;
twinN  = C.colorPre;

% Then we fit the preamble line to the real captures and to their twins.

[sR, iR, rR] = preambleLaw(sigma, realN);
[sT, iT, rT] = preambleLaw(sigma, twinN);

if any(isnan([sR iR rR sT iT rT]))
    error('preambleLaw returned NaN. Run testPreambleLaw to check it.');
end

% Next we print both fits and how far each capture sits off its line.

lineR   = sR * sigma + iR;
residR  = 100 * (realN ./ lineR - 1);
lineT   = sT * sigma + iT;
residT  = 100 * (twinN ./ lineT - 1);

fprintf('Real seven:   slope %9.0f /LSB  intercept %9.0f  r %.4f\n', sR, iR, rR);
fprintf('Colored twin: slope %9.0f /LSB  intercept %9.0f  r %.4f\n', sT, iT, rT);
fprintf('Twin slope / real slope = %.3f\n\n', sT / sR);
fprintf('%-14s %7s %9s %9s %8s %9s\n', 'capture', 'sigma', 'real', 'real fit', 'resid %', 'twin');

for k = 1:numel(sigma)
    fprintf('%-14s %7.4f %9d %9.0f %+8.2f %9d\n', A.capture(k), sigma(k), realN(k), lineR(k), ...
            residR(k), twinN(k));
end

fprintf('Worst real residual: %.2f %%\n\n', max(abs(residR)));

% After that we extend the real line to the twins with added noise to see how
% badly it misses.

tw = D(D.src == "twin", :);
lineD  = sR * tw.cellDC + iR;
residD = 100 * (tw.preambles ./ lineD - 1);

fprintf('Real line extrapolated to dithered twins:\n');
fprintf('%-14s %7s %9s %9s %8s\n', 'source', 'sigma', 'twin', 'line', 'twin/line');

for k = 1:height(tw)
    fprintf('%-14s %7.3f %9d %9.0f %8.3f\n', tw.capture(k), tw.cellDC(k), tw.preambles(k), ...
            lineD(k), tw.preambles(k) / lineD(k));
end

% Lastly, we save every point and both fits.

grp = [repmat("real", numel(sigma), 1); repmat("colorTwin", numel(sigma), 1); ...
       "ditherTwin " + tw.capture];

P = table(grp, [A.capture; C.capture; tw.capture], [sigma; sigma; tw.cellDC], ...
          [realN; twinN; tw.preambles], [lineR; lineR; lineD], [residR; 100 * (twinN ./ lineR - 1); residD], ...
          'VariableNames', {'set', 'capture', 'sigma', 'preambles', 'realLine', 'pctOffRealLine'});
writetable(P, fullfile(cfg.outDir, 'preamble_law.csv'));

Fit = table(["real"; "colorTwin"], [sR; sT], [iR; iT], [rR; rT], [max(abs(residR)); max(abs(residT))], ...
            'VariableNames', {'set', 'slope', 'intercept', 'r', 'worstResidPct'});
writetable(Fit, fullfile(cfg.outDir, 'preamble_law_fit.csv'));

fprintf('\nWrote out/preamble_law.csv and out/preamble_law_fit.csv\n');
