% First we load every table the figures need and make the folder they go in.

cfg = expConfig();
figDir = fullfile(cfg.outDir, 'fig');
if ~exist(figDir, 'dir'), mkdir(figDir); end

rd = @(f) readtable(fullfile(cfg.outDir, f), 'TextType', 'string');
A = rd('audit.csv');  F = rd('fix_compare.csv');  C = rd('color_twins.csv');
D = rd('dither.csv'); J = rd('injection.csv');    S = rd('injection_summary.csv');
O = rd('occupancy.csv');

% Then we set up the colors, short capture names, and this receiver's sigma
% range in LSB (one step of the 8-bit sampler).

ink  = [0.15 0.15 0.15];  blue = [0 0.447 0.698];  orng = [0.902 0.624 0];
verm = [0.835 0.369 0];   grn  = [0 0.620 0.451];  sky  = [0.337 0.706 0.914];
purp = [0.8 0.475 0.655]; grey = [0.6 0.6 0.6];

rigBand = [0.7096 0.8890];
piaware = 4.05;

short = extractBefore(A.capture, 4);
isCtl = contains(A.capture, "control");

% For figure 1, we plot the noise dump1090 reports next to the true background.

fig = newFig(760, 420);
ax = axes(fig);  hold(ax, 'on');  box(ax, 'on');

x = (1:height(A))';
bg  = A.bgModelDb;
rep = A.bgModelDb + A.inflationDb;

for k = 1:numel(x)
    plot(ax, [x(k) x(k)], [bg(k) rep(k)], '-', 'Color', grey, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(ax, x(k) + 0.12, rep(k) + 0.1, sprintf('%+.2f', A.inflationDb(k)), 'FontSize', 9, 'Color', ink);
end

hB = plot(ax, x, bg, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', blue, 'LineWidth', 1.5);
hR = plot(ax, x, rep, 'o', 'MarkerSize', 8, 'MarkerFaceColor', verm, 'MarkerEdgeColor', verm);

patch(ax, [4.6 7.6 7.6 4.6], [-42 -42 -38.8 -38.8], sky, 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
      'HandleVisibility', 'off');
text(ax, 4.7, -38.95, 'triad: 005 / 006 / 007, all within 2 min', 'FontSize', 9, 'Color', ink);

lbl = short;  lbl(isCtl) = lbl(isCtl) + "\newline(1082 MHz)";
set(ax, 'XTick', x, 'XTickLabel', lbl, 'TickLabelInterpreter', 'tex', 'XLim', [0.5 7.7], ...
    'YLim', [-42 -38.6], 'FontSize', 10);

ylabel(ax, 'dBFS (uncalibrated)');
xlabel(ax, 'capture (controls tuned off-channel, no ADS-B traffic)');
legend(ax, [hR hB], {'dump1090 reported noise', 'background from the innermost cell (cellDC)'}, ...
       'Location', 'southwest', 'Box', 'off');
title(ax, 'Reported noise includes undecoded traffic: +0.24 to +0.98 dB here', 'FontWeight', 'normal');

exportgraphics(fig, fullfile(figDir, 'fig1_audit.png'), 'Resolution', 200);  close(fig);

% For figure 2, we plot the preamble count against the background, zoomed in on
% the left and out to 4 LSB on the right.

haveLine = false;

try
    [sR, iR, rR] = preambleLaw(A.sigmaDC, F.preambles);
    haveLine = ~any(isnan([sR iR rR]));
catch
end

fig = newFig(980, 420);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

tw3 = D(D.src == "twin" & startsWith(D.capture, "003"), :);
tw6 = D(D.src == "twin" & startsWith(D.capture, "006"), :);
re3 = D(D.src == "real" & startsWith(D.capture, "003"), :);
re6 = D(D.src == "real" & startsWith(D.capture, "006"), :);

for panel = 1:2
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');

    if panel == 1
        xl = [0.69 0.91];

    else
        xl = [0 4.4];

        patch(ax, [rigBand fliplr(rigBand)], [0 0 900 900], sky, 'FaceAlpha', 0.15, ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
        xline(ax, piaware, ':', 'PiAware target', 'Color', ink, 'LabelVerticalAlignment', 'middle', ...
              'HandleVisibility', 'off');

        h3 = plot(ax, tw3.cellDC, tw3.preambles / 1e3, '-s', 'Color', orng, 'MarkerFaceColor', orng);
        h6 = plot(ax, tw6.cellDC, tw6.preambles / 1e3, '-d', 'Color', purp, 'MarkerFaceColor', purp);
        hd = plot(ax, [re3.cellDC; re6.cellDC], [re3.preambles; re6.preambles] / 1e3, 'x', ...
                  'Color', ink, 'MarkerSize', 9, 'LineWidth', 1.2);

    end

    if haveLine
        xs = linspace(xl(1), xl(2), 50);
        hL = plot(ax, xs, (sR * xs + iR) / 1e3, '--', 'Color', grey, 'LineWidth', 1.2);
    end

    hT = plot(ax, A.sigmaDC, C.colorPre / 1e3, 'o', 'MarkerEdgeColor', blue, 'LineWidth', 1.2);
    hRe = plot(ax, A.sigmaDC, F.preambles / 1e3, 'o', 'MarkerFaceColor', verm, 'MarkerEdgeColor', verm);

    if panel == 1
        dx = [0.004 0.004 0.004 0.004 0.004 -0.013 0.003]';  dy = [-6 -6 -6 -6 -2 4 -9]';
        text(ax, A.sigmaDC + dx, F.preambles / 1e3 + dy, short, 'FontSize', 8, 'Color', ink);

        ylim(ax, [320 520]);
        title(ax, 'This rig: 0.71 to 0.89 LSB', 'FontWeight', 'normal');

        leg = {hRe, hT};  names = {'real capture', 'colored Gaussian twin, same sigma'};

        if haveLine
            leg{end + 1} = hL;  names{end + 1} = sprintf('line fit to real, r = %.3f', rR);
        end

        legend(ax, [leg{:}], names, 'Location', 'northwest', 'Box', 'off');

    else
        ylim(ax, [0 900]);
        title(ax, 'Dithered to 4 LSB: the count saturates', 'FontWeight', 'normal');
        legend(ax, [h3 h6 hd], {'twin of 003, dithered', 'twin of 006, dithered', ...
               'real 003 and 006, dithered'}, 'Location', 'south', 'Box', 'off');

    end

    xlim(ax, xl);  set(ax, 'FontSize', 10);
    xlabel(ax, 'background sigma, LSB (cellDC)');

end

ylabel(tl, 'dump1090 preambles per 30 s (thousands)');
exportgraphics(fig, fullfile(figDir, 'fig2_preambles.png'), 'Resolution', 200);  close(fig);

% For figure 3, we plot how often messages decode against SNR (signal-to-noise
% ratio).

fig = newFig(980, 420);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');

show = [0.5 0.72 1 2 6];  cols = {verm, orng, grn, sky, blue};
hs = gobjects(1, numel(show) + 1);  names = cell(1, numel(show) + 1);

for k = 1:numel(show)
    j = J(J.source == "synthetic" & abs(J.sigma - show(k)) < 1e-9, :);
    hs(k) = plot(ax, j.snrDb, j.pDecode, '-', 'Color', cols{k}, 'LineWidth', 1.5);
    names{k} = sprintf('%.2g LSB', show(k));
end

j = J(J.source == "real 006", :);
hs(end) = plot(ax, j.snrDb, j.pDecode, 'k--', 'LineWidth', 1.2);
names{end} = 'real 006 (quantized twice)';

set(ax, 'FontSize', 10, 'XLim', [6 20], 'YLim', [0 1.02]);
xlabel(ax, 'SNR before quantization, dB  (A^2 / 2\sigma^2)');
ylabel(ax, 'fraction of DF17 frames decoded (--fix)');

legend(ax, hs, names, 'Location', 'southeast', 'Box', 'off');
title(ax, 'Same frames, same analog SNR, different background', 'FontWeight', 'normal');

% On the right, we plot the extra SNR needed compared with a background of 6
% LSB.

ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
s = S(S.source == "synthetic", :);

patch(ax, [rigBand fliplr(rigBand)], [-0.3 -0.3 2.8 2.8], sky, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
      'HandleVisibility', 'off');
text(ax, 0.73, 2.6, 'this rig', 'FontSize', 9, 'Color', ink);
xline(ax, piaware, ':', 'PiAware target', 'Color', ink, 'LabelVerticalAlignment', 'top', ...
      'HandleVisibility', 'off');
yline(ax, 0, '-', 'Color', grey, 'HandleVisibility', 'off');

h50 = plot(ax, s.sigma, s.loss50, '-o', 'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.5);
h90 = plot(ax, s.sigma, s.loss90, '-s', 'Color', verm, 'MarkerFaceColor', verm, 'LineWidth', 1.5);

% Then we add a second row of tick labels: the noise dump1090 would report with
% no traffic.

ticks = [0.5 0.72 1 1.5 2 3 4 6];
lut = dumpLUT();  mu = 127.397;  tl2 = strings(size(ticks));

for k = 1:numel(ticks)
    p = quantModel(mu, ticks(k));
    tl2(k) = sprintf('%.4g\\newline%.1f', ticks(k), 10 * log10(sum(sum((p * p') .* lut.^2)) / 65535^2));
end

set(ax, 'XScale', 'log', 'XTick', ticks, 'XTickLabel', tl2, 'TickLabelInterpreter', 'tex', ...
    'XLim', [0.45 6.6], 'YLim', [-0.3 2.8], 'FontSize', 10);
xlabel(ax, {'background sigma, LSB', 'reported noise with no traffic, dBFS'});
ylabel(ax, 'extra SNR needed vs 6 LSB, dB');

legend(ax, [h50 h90], {'at 50 % decode', 'at 90 % decode'}, 'Location', 'north', 'Box', 'off');
title(ax, 'Quantization loss (synthetic, 900 frames per 1 dB)', 'FontWeight', 'normal');

exportgraphics(fig, fullfile(figDir, 'fig3_injection.png'), 'Resolution', 200);  close(fig);

% For figure 4, we plot how much each estimate rises as messages cover more of
% the file.

fig = newFig(980, 420);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for sg = [0.72 4]
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    o = sortrows(O(abs(O.sigma - sg) < 1e-9, :), 'occPct');

    h1 = plot(ax, o.occPct, o.dumpExcessDb, '-o', 'Color', verm, 'MarkerFaceColor', verm, 'LineWidth', 1.5);
    h2 = plot(ax, o.occPct, o.cellDCDb, '-s', 'Color', orng, 'MarkerFaceColor', orng, 'LineWidth', 1.5);
    h3 = plot(ax, o.occPct, o.centralFitDb, '-d', 'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.5);
    h4 = plot(ax, o.occPct, o.centralFitOracleDb, ':', 'Color', grey, 'LineWidth', 1.5);

    set(ax, 'FontSize', 10, 'XLim', [0 20.5], 'YLim', [-0.1 2.4]);
    xlabel(ax, 'samples covered by injected frames, %');
    title(ax, sprintf('background %.2g LSB', sg), 'FontWeight', 'normal');

    if sg == 0.72
        ylabel(ax, 'rise over the same file at 0 %, dB');
        legend(ax, [h1 h2 h3 h4], {'dump1090 reported noise', 'cellDC', 'centralFit', ...
               'centralFit, traffic samples removed'}, 'Location', 'northwest', 'Box', 'off');
    end
end

title(tl, 'Synthetic busy site: both candidates drift with traffic, less than the reported noise', 'FontSize', 11);
exportgraphics(fig, fullfile(figDir, 'fig4_occupancy.png'), 'Resolution', 200);  close(fig);

% Lastly, we print where the figures went.

fprintf('Wrote 4 figures to %s (preamble line drawn: %d)\n', figDir, haveLine);

function fig = newFig(w, h)

% Here we make a hidden figure with a light background.

fig = figure('Visible', 'off', 'Position', [0 0 w h]);

try
    theme(fig, 'light');
catch
    set(fig, 'Color', 'w');
end

end
