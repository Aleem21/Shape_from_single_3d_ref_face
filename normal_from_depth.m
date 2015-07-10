function [ n,N_ref ] = normal_from_depth( depth_map )
%NORMAL_FROM_DEPTH Summary of this function goes here
%   Detailed explanation goes here
delta = [2 2];

% order of 1,-1 is flipped because conv2 function adds another flip itself.
fx = [0 -1 1];
fy = [0 -1 1]';
fx = [-1 0 1];
fy = [-1 0 1]';
p = conv2(depth_map, fx, 'same')/delta(1);
q = conv2(depth_map, fy, 'same')/delta(2);
N_ref = (p.^2 + q.^2 +1).^0.5;
n(:,:,1) = p./N_ref;
n(:,:,2) = q./N_ref;
n(:,:,3) = -1./N_ref;


n(end,:,:) = NaN;
n(:,end,:) = NaN;
N_ref(end,:) = NaN;
N_ref(:,end) =NaN;

n(1,:,:) = NaN;
n(:,1,:) = NaN;
N_ref(1,:) = NaN;
N_ref(:,1) =NaN;

end

