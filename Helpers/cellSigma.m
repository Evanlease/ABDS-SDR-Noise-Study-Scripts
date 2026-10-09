function sigma = cellSigma(P)

% Here we turn the share of samples in the 4 center squares into sigma, the
% spread of the noise.

sigma = 1 ./ (sqrt(2) * erfinv(sqrt(P)));

end
