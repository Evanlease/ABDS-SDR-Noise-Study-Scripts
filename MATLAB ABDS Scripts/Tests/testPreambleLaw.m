% First we set up the sigmas, the preamble counts, and the answers we expect.

sigma     = [0.8486; 0.7665; 0.7096; 0.8890; 0.7207; 0.7174; 0.7177];
preambles = [461185; 400498; 342520; 498369; 357695; 359962; 355392];
expect    = struct('slope', 829519.4, 'intercept', -239783.4, 'r', 0.998);

% Then we run preambleLaw and check the slope, intercept and r.

[slope, intercept, r] = preambleLaw(sigma, preambles);
ok = [abs(slope - expect.slope) < 1, abs(intercept - expect.intercept) < 1, ...
      abs(r - expect.r) < 5e-4];

fprintf('slope     %12.1f  expected %12.1f  %s\n', slope, expect.slope, passFail(ok(1)));
fprintf('intercept %12.1f  expected %12.1f  %s\n', intercept, expect.intercept, passFail(ok(2)));
fprintf('r         %12.4f  expected %12.3f  %s\n', r, expect.r, passFail(ok(3)));

% Lastly, we start an if statement that stops with an error if any check failed.

if ~all(ok)
    error('preambleLaw: FAIL, see above.');
end

fprintf('\nPASS. Now run runPreambleLaw.\n');

function s = passFail(ok)

% Here we turn true or false into PASS or FAIL.

if ok, s = 'PASS'; else, s = 'FAIL'; end

end
