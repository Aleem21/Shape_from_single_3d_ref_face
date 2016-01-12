function [ out ] = scale_up_nan( inp,sz,nans )
%SCALE_DOWN_NAN Summary of this function goes here
%   Detailed explanation goes here

small = inp;
small(isnan(small)) = 0;

nan_small = double(~(isnan(inp)));
big = impyramid(small,'expand');
% nan_small = imresize(nan_big,scale);
nan_big = impyramid(nan_small,'expand');
out = big./nan_big;
out(isnan(out)) = 0;
out = imresize(out,sz);
out(nans) = nan;

end

