function [ l ] = fit_lighting( mData)
%FIT_LIGHTING Summary of this function goes here
%   Detailed explanation goes here
alb = mData(1,:);
im = mData(2,:);
Y = mData(3:end,:);

Y_alb = Y.*repmat(alb,size(Y,1),1);
l = im * pinv(Y_alb);

end

