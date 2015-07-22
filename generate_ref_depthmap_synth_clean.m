function [ depth_map, n, N_ref,pts ] = generate_ref_depthmap_synth_clean( rrange,sigma,scale,talk )
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
if nargin < 3
    scale = 1;
end
if nargin < 4
    talk = 0;
end

import plot3D_helper.plot_pointcloud
import plot3D_helper.label_axis
import plyread.plyread
import hidden_point_removal.HPR
[r,~] = meshgrid(0:rrange,0:rrange);
gauss = fspecial('gaussian',size(r),sigma);
gauss = gauss/max(gauss(:));
depth_map = gauss/max(gauss(:))*5*scale;
depth_map = max(depth_map(:))+1-depth_map;

depth_map(:,1) = depth_map(:,2);
depth_map(:,end-1) = depth_map(:,end-2);
depth_map(:,end) = depth_map(:,end-2);


%% find normals from depth
[ n,N_ref ] = normal_from_depth( depth_map );

depth_map(:,end) = NaN;
depth_map(end,:) = NaN;
depth_map(:,1) = NaN;
depth_map(1,:) = NaN;




if talk
    figure;
    subplot(3,3,[1,2, 4,5 , 7,8])
    imagesc(depth_map)
    xlabel('y')
    ylabel('x')
    title('Depth map')
    subplot(3,3,3)
    imagesc(n(:,:,1))
    title('n_x')
    subplot(3,3,6)
    imagesc(n(:,:,2))
    title('n_y')
    subplot(3,3,9)
    imagesc(n(:,:,3))
    title('n_z')
    truesize
end
