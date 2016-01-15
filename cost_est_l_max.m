function [ out ] = cost_est_l_max( l,radiosity,Y )
%COST_EST_L Summary of this function goes here
%   Detailed explanation goes here
% out = sum((l'*Y - radiosity).^2);
l = l';
out = max(l(1)*Y(1,:) + l(2:4)*Y(2:4,:),0) - radiosity;
% out = l(1)*Y(1,:) + l(2:9)*Y(2:9,:) - radiosity;

out(isinf(out))=0;
out(isnan(out))=0;
end

