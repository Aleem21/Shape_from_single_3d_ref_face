function [ boundary,in_bound ] = find_boundary( bw,full )
%FIND_BOUNDARY Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    full = 1;
end
if full
    box = ones(3);
else
    box = [0 1 0;
        1 0 1;
        0 1 0];
end
convbw = conv2(double(bw),box,'same');
boundNout = convbw<sum(box(:));
boundary = boundNout.*bw;
if nargout >1
    inside = (bw-boundary)>0;
    in_bound = find_boundary(inside,full);
end
end

