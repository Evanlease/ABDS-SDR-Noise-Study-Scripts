function [slope, intercept, r] = preambleLaw(sigma, preambles)

% First we fit a straight line to the preamble counts against sigma.

c = polyfit(sigma, preambles, 1);
slope = c(1);
intercept = c(2);

% Lastly, we find r, which says how close the points are to a straight line.

R = corrcoef(sigma, preambles);
r = R(1, 2);

end
