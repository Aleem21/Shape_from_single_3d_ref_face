function [ cost ] = cost_lighting( l,mData )
%COST_LIGHTING Summary of this function goes here
%   Detailed explanation goes here
alb = mData(1,:);
im = mData(2,:);
Y = mData(3:end,:);

cost = max(l*Y.*alb,0) - im;
end

