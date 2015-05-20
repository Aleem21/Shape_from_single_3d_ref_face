function [ out ] = find_elements( list,r,c )
%FIND_ELEMENTS Summary of this function goes here
%   Detailed explanation goes here
[rs] = find(list(:,1)==r);
cs_ind = list(rs,2)==c;

out = rs(cs_ind);

end

