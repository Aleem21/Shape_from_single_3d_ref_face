function [ mew,sigma ] = compute_gaussian_mixture( x,regions )
%COMPUTE_GAUSSIAN_MIXTURE Summary of this function goes here
%   Detailed explanation goes here
if size(x,1)>size(x,2)
    x = x';
end
r_min = min(regions);
r_max = max(regions);
r = r_min:r_max;
for i=1:numel(r)
    [mew{i}, sigma{i}] = gaussian(x(:,regions==r(i)));
end
end

