function [ n,N_ref, depth_map ] = normal_from_depth( depth_map )
%NORMAL_FROM_DEPTH Summary of this function goes here
%   Detailed explanation goes here
delta = [1 1];
fx = [1 -1]';
fy = [1 -1];

p = conv2(depth_map, fx, 'same')/delta(1);
q = conv2(depth_map, fy, 'same')/delta(2);
N = 1./(p.^2 + q.^2 +1).^0.5;
n(:,:,1) = N.*p;
n(:,:,2) = N.*q;
n(:,:,3) = N;
N_ref = 1./N;


depth_map(:,end) = NaN;
depth_map(end,:) = NaN;
n(end,:,:) = NaN;
n(:,end,:) = NaN;
N_ref(end,:) = NaN;
N_ref(:,end) =NaN;

end

