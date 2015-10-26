function [ n,N_ref ] = normal_from_depth( depth_map )
%NORMAL_FROM_DEPTH Summary of this function goes here
%   Detailed explanation goes here

% order of 1,-1 is flipped because conv2 function adds another flip itself.
fx= [1 -1];
fy = [1 -1]';
delta = [1 1];

depth_map2 = depth_map(end:-1:1,end:-1:1);
p = conv2(depth_map2, fx, 'same')/delta(1);
q = conv2(depth_map2, fy, 'same')/delta(2);
p = p(end:-1:1,end:-1:1);
q = q(end:-1:1,end:-1:1);

% 
% fx = [0 -1 1];
% fy = [0 -1 1]';
% delta = [1 1];
% % fx = [-1 0 1];
% % fy = [-1 0 1]';
% % delta = [2 2];
% 
% p = conv2(depth_map, fx, 'same')/delta(1);
% q = conv2(depth_map, fy, 'same')/delta(2);

p2 = p(:,2:end);
p2(:,end+1)=0;
p(:,1) = p2(:,1);
p(isnan(p)) = p2(isnan(p));

q2 = q(2:end,:);
q2(end+1,:)=0;
q(1,:) = q2(1,:);
q(isnan(q)) = q2(isnan(q));


N_ref = (p.^2 + q.^2 +1).^0.5;
n(:,:,1) = p./N_ref;
n(:,:,2) = q./N_ref;
n(:,:,3) = -1./N_ref;

% 
% n(end,:,:) = NaN;
% n(:,end,:) = NaN;
% N_ref(end,:) = NaN;
% N_ref(:,end) =NaN;
% 
% n(1,:,:) = NaN;
% n(:,1,:) = NaN;
% N_ref(1,:) = NaN;
% N_ref(:,1) =NaN;

end

