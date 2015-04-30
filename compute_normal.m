function [ n ] = compute_normal( depth_map )
%COMPUTE_NORMAL computes normal from depth_map. The output n is a three
%slice matrix with row and col dimentions same as depth_map

%% find normal
fx = [-1 0 1]';         %x is vertical, starting from top
fy = [-1 0 1];          %y is horizontal, starting from left

p = conv2(depth_map, fx, 'same');
q = conv2(depth_map, fy, 'same');
weight = 1./(p.^2 + q.^2 +1).^0.5;
n(:,:,1) = weight.*p;
n(:,:,2) = weight.*q;
n(:,:,3) = -weight;

end

