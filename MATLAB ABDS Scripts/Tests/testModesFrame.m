% First we build a known example message with modesFrame.

[~, hex] = modesFrame(hex2dec('4840D6'), 'KLM1023');

% Lastly, we check that it matches the published answer and print the result.

ok = hex == "8D4840D6202CC371C32CE0576098";
fprintf('modesFrame reproduces 8D4840D6202CC371C32CE0576098: %s (got %s)\n', ...
        string(ok), hex);
