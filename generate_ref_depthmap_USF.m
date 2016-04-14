function [ depth_map, n, N_ref,albedo, eye_map,scalez,pts,model ] = generate_ref_depthmap_USF( Scale, Rpose, im,im_c, dataset_path,use_model,talk )

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
RGB_to_XYZ = [  0.4124564 0.3575761 0.1804375;
    0.2126729 0.7151522 0.0721750;
    0.0193339 0.1191920 0.9503041];
XYZ_to_RGB = inv(RGB_to_XYZ);
import plot3D_helper.plot_pointcloud
import plot3D_helper.label_axis
import plyhelper.plyread
if nargin<4
    use_model = 0;
end
if nargin < 5
    talk = 0;
end


xrange = [1 size(im,2)];
yrange = [1 size(im,1)]; 

%% Generate reference pointcloud and texture map (albedo)
black_eyes = 1;
is_full = 0;
[pts,tri,rgb,~,~,~,~,~,eye_rgb] = read_USF_eko(dataset_path,512,512,black_eyes,talk,is_full);
%% data conditioning

%pose correction
pts_rotated = Rpose*[pts([1 2 3],:); ones(1,size(pts,2))];
pts_rotated = pts_rotated(1:3,:);
pts_rotated(3,:) = pts(3,:)*min(Scale);
pts_rotated = pts_rotated([1 2 3],:);


%% generate depth map by calling mex file
if use_model
    [FV,b,R,t,s] = fit_model(im_c,1);
    model.FV = FV; model.b = b; model.R = R; model.t = t; model.s = s;
    % get depth
    [pts_rotated,tri,clr] = model2GL(model,mean(xrange));
    model.clr = clr;
    [depth_map,~] = computer_depth_USF( pts_rotated,tri,clr,xrange,yrange,im_c,talk );
    % fill in the gap between lips
    steps = abs(conv2(depth_map,[-1 1]','same'))>45;
    steps([1 end],:) = 0;
    steps(:,[1 end]) = 0;
    [r,c] = find(steps);
    for i = 1:numel(r)
        depth_map(r(i),c(i)) = median(depth_map(sub2ind(size(depth_map),r(i)+[-1 0 1],c(i)+[0 0 0])));
    end
    
    %get albedo
    load('./Basel Face Model/01_MorphableModel.mat');
    clr = double(reshape(texMU,3,[]))/255;
	[~,albedo] = computer_depth_USF( pts_rotated,tri,clr,xrange,yrange,im_c,talk );
    
else
    [depth_map,albedo] = computer_depth_USF( pts_rotated,tri,rgb,xrange,yrange,im_c,talk );
end
scalez = 1/9*1/min(Scale);
% albedo = correct(albedo,RGB_to_XYZ);
% albedo = double(albedo(:,:,2));
% albedo = correct(albedo,RGB_to_XYZ);
% albedo = xyz2rgb(albedo);
        albedo = double(rgb2gray(albedo));
% albedo = double(albedo(:,:,2));
albedo = albedo/max(albedo(:));
%% generate eyemap
if use_model
    load('eye_mask.mat');
    clr = double(repmat(eye_mask',3,1)>250);
	[~,eye_map] = computer_depth_USF( pts_rotated,tri,clr,xrange,yrange,im_c,0 );
    eye_map= eye_map(:,:,1)~=1;
else
    [~,eye_map] = computer_depth_USF( pts_rotated,tri,eye_rgb,xrange,yrange,im_c,talk );
    eye_map = double(eye_map(:,:,1)~=1);
end
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
