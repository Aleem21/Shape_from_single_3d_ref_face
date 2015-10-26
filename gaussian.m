function [ mew,sigma ] = gaussian( x )
%COMPUTE_GAUSSIAN Summary of this function goes here
%   Detailed explanation goes here
sz = size(x);
if sz(1)>sz(2)
x = x';
end
mew = mean(x,2);
centered = x-repmat(mew,1,size(x,2));
sigma = (centered*centered')/size(x,2);
