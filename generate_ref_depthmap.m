function [ depth_map, n, N_ref,pts ] = generate_ref_depthmap( ply_path,Scale, talk,rRes,cRes,Rpose,im,center,xrange,yrange )
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
if nargin < 6
    Rpose = eye(4);
end
if nargin < 7
    im = [];
end
if nargin < 9
    xrange = [-1 1];
end
if nargin < 10
    yrange = [-1 1];
end
if nargin <8
    center = 1;
end
%% read model
[tri,pts] = plyread(ply_path,'tri');

%% data conditioning
pts = [pts ones(size(pts,1),1)]';


% Rpose = makehgtform('scale',0.005);
%pose correction
pts = Rpose*pts;
pts = pts(1:3,:);

%% Draw point cloud
if talk > 2
    figure;
    plot_pointcloud(pts,'.b')
    label_axis();
end

%% generate depth map by calling mex file
depth_map = compute_depthmap( pts, tri, xrange, yrange, cRes, rRes,center )*min(Scale);

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
% width of filter in x and y direction
% delta = [2 2];
% fx = [-1 0 1]';         %x is vertical, starting from top
% fy = [-1 0 1];          %y is horizontal, starting from left

delta = [1 1];
fx = [-1 1]';
fy = [-1 1];

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

if talk
    figure;
    subplot(3,3,[1,2, 4,5 , 7,8])
    imshow(n/2+0.5)
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
