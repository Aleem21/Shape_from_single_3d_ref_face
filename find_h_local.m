function [ N_H_theta,H ] = find_h_local( N,w_o,wi_grid )
%FIND_H_LOCAL Finds half angles for wi_grid, when wi_grid is expressed in
%local frame. The local frame is such that its [0;0;1] appears to be
%pointing in N direction in global frame. w_o is defined in global
%frame.
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

%% Contruct rotation matrix to rotate w_o into local frame
N = N/norm(N);
R(3,:) = N;
temp = cross(R(3,:)',rand(3,1));
R(2,:) = temp/norm(temp);
R(1,:) = cross(R(2,:),R(3,:));
if det(R)<0
    R(1:2,:) = R([2 1],:);
end
% Do the rotation
w_o_local = R*w_o;
% Find H in local coordinates
grid_sph = aziElev2aziIncl(wi_grid);
[x,y,z] = sph2cart(grid_sph(:,1),grid_sph(:,2),ones(size(grid_sph,1),1));
grid_cart = [x';y';z'];
H = grid_cart + repmat(w_o_local,1,size(grid_cart,2));
H = H./repmat(sqrt(sum(H.^2)),3,1);
N_H_theta = acos([0;0;1]'*H)';
end

