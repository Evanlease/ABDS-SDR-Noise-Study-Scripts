function [bI, bQ] = synthNoise(sigma, N, seed, mu, qRatio)

% First we use this receiver's center and equal I and Q if they weren't given.

if nargin < 4 || isempty(mu), mu = 127.397; end
if nargin < 5 || isempty(qRatio), qRatio = 1; end

% Then we create random noise, round it to bytes, and keep it between 0 and 255.

rng(seed, 'twister');
bI = uint8(min(max(round(mu + sigma * randn(N, 1)), 0), 255));
bQ = uint8(min(max(round(mu + qRatio * sigma * randn(N, 1)), 0), 255));

end
