function [ depth_map, n, N_ref,pts ] = generate_ref_depthmap( ply_path, rRes, cRes, Rpose, im, xrange, yrange, talk )
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
%       talk = 1 means only plot depthmap imagesc alignment and normals
%       talk = 2 means also plot depthmap surf
%       talk = 3 means also draw input ref model
%
%   Written by Muhammad Ahmed Riaz(mriaz@ucsd.edu) on 15th April, 2015.

import plot3D_helper.plot_pointcloud
import plot3D_helper.label_axis
import plyread.plyread
import hidden_point_removal.HPR
if nargin <8
    talk = 0;
end
if nargin < 7
    yrange = [-1 1];
end
if nargin < 6
    xrange = [-1 1];
end
if nargin < 5
    im = [];
end
if nargin < 4
    Rpose = eye(4);
end
if nargin < 3
    cRes = 400;
end
if nargin < 2
    rRes = 400;
end

%% read model
[tri,pts] = plyread(ply_path,'tri');

%% data conditioning
pts = [pts ones(size(pts,1),1)]';

%pose correction
pts = Rpose*pts;
pts = pts(1:3,:);



%% generate depth map by calling mex file
depth_map = compute_general_depthmap( pts, tri, cRes, rRes, xrange, yrange);
% depth_map = compute_depthmap( pts, tri, xrange, yrange, cRes, rRes,center )*min(Scale);

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

%% find normals from depth
[ n,N_ref ] = normal_from_depth( depth_map );

depth_map(:,end) = NaN;
depth_map(end,:) = NaN;


%% draw stuff

% Draw point cloud
if talk > 2
    figure;
    plot_pointcloud(pts,'.b')
    label_axis();
end

if talk
    figure;
    subplot(3,3,[1,2, 4,5 , 7,8])
    imagesc(depth_map)
    xlabel('x')
    ylabel('y')
    title('depth map')
    axis equal; axis off
    subplot(3,3,3)
    imagesc(n(:,:,1))
    title('n_x')
    axis equal; axis off
    subplot(3,3,6)
    imagesc(n(:,:,2))
    title('n_y')
    axis equal; axis off
    subplot(3,3,9)
    imagesc(n(:,:,3))
    title('n_z')
    axis equal; axis off
    truesize
end
