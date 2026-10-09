function H = byteHist(bI, bQ)

% Here we count how many samples land on each pair of I and Q bytes.

H = accumarray([double(bI) + 1, double(bQ) + 1], 1, [256 256]);

end
