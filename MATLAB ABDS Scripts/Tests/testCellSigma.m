% First we set up the two controls' center shares and the sigmas we expect.

P       = [0.6928; 0.5420];
expect  = [0.7248; 0.8949];
cellDC  = [0.7174; 0.8890];
names   = {'006 control'; '004 control'};

% Then we run cellSigma and check each answer.

got = cellSigma(P);
ok = abs(got - expect) < 5e-4;

for k = 1:numel(P)
    fprintf('%-12s P = %.4f  sigma = %.4f  expected %.4f  %s\n', ...
            names{k}, P(k), got(k), expect(k), passFail(ok(k)));
end

% Lastly, we start an if statement where, if both pass, we print how far each
% answer is from cellDC in dB.

if all(ok)
    fprintf('\nPASS. Centred minus cellDC: %+.3f dB (006), %+.3f dB (004).\n', ...
            20 * log10(got ./ cellDC));
    fprintf('That is what assuming the centre is 127.5 instead of 127.397 costs.\n');
else
    error('cellSigma: FAIL, see above.');
end

function s = passFail(ok)

% Here we turn true or false into PASS or FAIL.

if ok, s = 'PASS'; else, s = 'FAIL'; end

end
