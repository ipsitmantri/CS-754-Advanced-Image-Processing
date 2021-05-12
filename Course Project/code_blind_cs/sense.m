% function to sense the signal. Here y represents the sensed signal and 
% cnt represents a matrix which counts the number of times a pixel is
% measured and will be used for avergaing later
function [y, cnt] = sense(x, seed, ratio)
    [rows cols] = size(x);
    y = zeros(floor(ratio*64), (rows-7)*(cols-7));
    iter = 1;
    cnt = zeros(size(x));
    rng(seed);
    for i=1:rows-7
        for j=1:cols-7
            phi = randn(floor(ratio*64), 64);
            patch = x(i:i+7, j:j+7);
            y(:, iter) = phi*patch(:);
            iter=iter+1;
            cnt(i:i+7, j:j+7) = cnt(i:i+7, j:j+7)+1;
        end
    end
end