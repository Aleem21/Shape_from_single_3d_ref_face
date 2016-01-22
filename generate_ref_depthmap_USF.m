function [ depth_map, n, N_ref,albedo, eye_map,scalez,pts ] = generate_ref_depthmap_USF( Scale, Rpose, im,im_c, dataset_path,talk )

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

if nargin < 4
    talk = 0;
end


xrange = [1 size(im,2)];
yrange = [1 size(im,1)]; 

%% Generate reference pointcloud and texture map (albedo)
black_eyes = 1;
[pts,tri,rgb,~,~,~,~,~,eye_rgb] = read_USF_eko(dataset_path,512,512,black_eyes,talk);
%% data conditioning

%pose correction
pts_rotated = Rpose*[pts([1 2 3],:); ones(1,size(pts,2))];
pts_rotated = pts_rotated(1:3,:);
pts_rotated(3,:) = pts(3,:)*min(Scale);
pts_rotated = pts_rotated([1 2 3],:);


%% generate depth map by calling mex file
[depth_map,albedo] = computer_depth_USF( pts_rotated,tri,rgb,xrange,yrange,im_c,talk );
scalez = 1/9*1/min(Scale);
albedo = rgb2xyz(albedo);
albedo = double(albedo(:,:,2));

%% generate eyemap
[~,eye_map] = computer_depth_USF( pts_rotated,tri,eye_rgb,xrange,yrange,im_c,talk );
eye_map = double(eye_map(:,:,1)~=1);
% eye_map(eye_map==0) = 0.1;
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
    plot(pts(1,:),size(im,1)-pts(2,:),'.r')
end

%% find normals from depth
[ n,N_ref ] = normal_from_depth( depth_map );

depth_map(:,end) = NaN;
depth_map(end,:) = NaN;
depth_map(:,1) = NaN;
depth_map(1,:) = NaN;

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
