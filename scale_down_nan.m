function [ out ] = scale_down_nan( inp )
%SCALE_DOWN_NAN Summary of this function goes here
%   Detailed explanation goes here

big = inp;
big(isnan(big)) = 0;

nan_big = double(~(isnan(inp)));
small = impyramid(big,'reduce');
% nan_small = imresize(nan_big,scale);
nan_small = impyramid(nan_big,'reduce');
out = small./nan_small;
end

