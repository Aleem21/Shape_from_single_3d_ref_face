function [ out ] = cost_est_l( l,radiosity,Y )
%COST_EST_L Summary of this function goes here
%   Detailed explanation goes here
% out = sum((l'*Y - radiosity).^2);
l = l';
out = sum((l(1)*Y(1,:) + max(l(2:4)*Y(2:4,:),0) - radiosity).^2)+size(Y,2)*l(1)/100;

end

