function [dI, dQ] = ditherIQ(bI, bQ, sigmaD, seed)

% Here we add random noise of size sigmaD to every byte, then round it and keep
% it between 0 and 255.

rng(seed, 'twister');
n = numel(bI);

dI = uint8(min(max(round(double(bI) + sigmaD * randn(n, 1)), 0), 255));
dQ = uint8(min(max(round(double(bQ) + sigmaD * randn(n, 1)), 0), 255));

end
