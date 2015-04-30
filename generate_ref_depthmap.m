function [ depth_map, n, pts ] = generate_ref_depthmap( ply_path,Scale, talk,rRes,cRes,Rpose,im,xrange,yrange )
%GENERATE_REF_DEPTHMAP generated depthmap from reading the ply file of the
%model, resolution of rows and columns specified by inputs
%
%   ply_path is the path of the ply file to read from
%
%   rRes,cRes   resolution of rows and columns in the depth map
%
%   xrange,yrange are of the form [minXorY maxXorY], which is the range in
%               the model for which the depth map is to be generated
%
%   talk is verbosity;
%       talk = 0 means no display
%       talk = 1 means only plot depthmap imagesc and alignment
%       talk = 2 means also plot depthmap surf
%       talk = 3 means also draw input ref model
%
%   Written by Muhammad Ahmed Riaz(mriaz@ucsd.edu) on 15th April, 2015.

import plot3D_helper.plot_pointcloud
import plot3D_helper.label_axis
import plyread.plyread
import hidden_point_removal.HPR
if nargin < 2
    Scale = [1 1];
end
if nargin <3
    talk = 0;
end
if nargin < 4
    rRes = 200;
end
if nargin < 5
    cRes = 200;
end
if nargin < 8
    xrange = [-1 1];
end
if nargin < 9
    yrange = [-1 1];
end

%% read model
[tri,pts] = plyread(ply_path,'tri');

%% data conditioning
pts = [pts ones(size(pts,1),1)]';

%pose correction
pts = Rpose * pts;
pts = pts(1:3,:);

%% Hidden Points Removal
% pts_old = pts;
% pts = pts_old;
% inds = hidden_point_removal.HPR(pts',[0 0 1], 1);
% pts = pts(:,inds);

%% Draw point cloud
if talk > 2
    figure;
    plot_pointcloud(pts,'.b')
    label_axis();
end

%% generate depth map by calling mex file
depth_map = compute_depthmap( pts, tri, xrange, yrange, cRes, rRes )*min(Scale);

%% draw depth map
if talk >0
    figure; imagesc(depth_map)
    label_axis(2)
end
if talk >1
    figure;
    surf(depth_map(end:-1:1,:),'edgealpha',0);
    label_axis(3)
    
    figure; imshow(im); hold on;
    plot(pts(1,:),rRes-pts(2,:),'.r')
end


%% find normal
%pix_width = 2./[rRes; cRes].*Scale(:);
%pix_width = 2./[rRes; cRes];
%delta = 2*pix_width;
% width of filter in x and y direction
delta = [2 2];
fx = [1 0 -1]';         %x is vertical, starting from top
fy = [1 0 -1];          %y is horizontal, starting from left

p = conv2(depth_map, fx, 'same')/delta(1);
q = conv2(depth_map, fy, 'same')/delta(2);
weight = 1./(p.^2 + q.^2 +1).^0.5;
n(:,:,1) = weight.*p;
n(:,:,2) = weight.*q;
n(:,:,3) = weight;
