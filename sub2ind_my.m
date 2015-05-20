function [ ind ] = sub2ind_my( sz, r,c )
%SUB2IND_MY Summary of this function goes here
%   Detailed explanation goes here

ind = r(:) + (c(:)-1)*sz(1);
end

