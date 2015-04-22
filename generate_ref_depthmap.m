function [ depth_map ] = generate_ref_depthmap( ply_path,talk,rRes,cRes,xrange,yrange )
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
%       talk = 1 means only plot depthmap
%       talk = 2 means draw both depthmap and the model being input
%
%   Written by Muhammad Ahmed Riaz(mriaz@ucsd.edu) on 15th April, 2015.

import plot3D_helper.plot_pointcloud
import plot3D_helper.label_axis
import plyread.plyread

if nargin <2
    talk = 0;
end
if nargin < 3
    rRes = 200;
end
if nargin < 4
    cRes = 200;
end
if nargin < 5
    xrange = [-1 1];
end
if nargin < 6
    yrange = [-1 1];
end

%% read model
[~,pts] = plyread(ply_path,'tri');

%% model specific conditioning
Rx = makehgtform('xrotate',deg2rad(20)); Rx = Rx(1:3,1:3);
Ry = makehgtform('yrotate',deg2rad(07)); Ry = Ry(1:3,1:3);
Rz = makehgtform('zrotate',deg2rad(-03)); Rz = Rz(1:3,1:3);

pts = Rz*Ry*Rx*pts';
pts(2,:) = pts(2,:)+0.35;

%% more general conditioning
pts = pts';

%% draw model point cloud
if talk > 1
    plot_pointcloud(pts,'.b')
    label_axis();
end

%% generate point cloud
[Xi,Yi] = meshgrid(linspace(xrange(1),xrange(2),cRes),linspace(yrange(1),yrange(2),rRes));
depth_map = griddata (pts(:,1),pts(:,2),pts(:,3),Xi,Yi);

%% draw depth map
if talk >0
    figure; surf(depth_map,'edgealpha',0)
    label_axis(2)
end

end

