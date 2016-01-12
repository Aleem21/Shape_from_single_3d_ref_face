clear variables
close all
import plot3D_helper.label_axis
fusion_path = '.\data\fusion.jpg';
addpath ('.\SHT')
%% set initial variables
% optimization weights
lambda1 = 50;
lambda2 = 30;

% number of maximum iterations of depth optimization
max_iter = 50;

% do you want to put in albedo into input image synthesis (keep it 1,
% option was created for early experiments to create the problems easier
% with no albedo complications
is_albedo = 1;

% Keep fixed. Legacy options
lambda_bound = 1;
is_dz_depth = 0;
lambda_dz = 0;

% Do you want to optimize for albedo as well?
is_alb_opt = 1;

% Use jacobian in cost function optimization for depth
jack = 'on';

% Compute both albedo and depth at the same time in a joint optimization?
% Keep 'off'. Option available because of old experiments, now obselete
combined = 0;

is_alb_dz = 1;

% Choose one of the two algorithms for optimizations. LM is faster than
% trust region but at times gives worse results
algo = 'levenberg-marquardt';
algo = 'trust-region-reflective';

% 1=old boundary condition. 2=new boundary condition. 2 is the better
% option and keep it to 2 for better results. option 1 not fully supported
% by all functions. Option exists for legacy reasons.
boundary_type = 2;

% To compute flow or not? Keep to 1 for better results
flow = 1;

% Use some real image as input rather than synthesizing one from dataset?
internet = 1;
internet_path = '.\data\internet\';
internet_imgs = {'paper small.png'};

raw = 1;
raw_path = '.\data\RAW\';
raw_imgs = {'mahmood2.cr2'};

% USF images path
folder_path = '.\data\USF_images\';

% change this for different USF images
impaths = {'03721c15.eko'};
n = numel(impaths);
f=figure;hold on
count = 1;

% talk = 4 - fully verbose
% talk = 0 - least talking and result displaying.
% talk is not fully consistent in all the code. Needs improvement.
talk = 4;

i=1;
impath = [folder_path impaths{i}];
%% make image
if raw
    im_c = imread_raw([raw_path raw_imgs{i}]);
    mask = double(imread([raw_path raw_imgs{i}(1:end-4) '_mask.tiff'])>0);
    mask = mask(:,:,1);
    im_c = imresize(im_c,0.07);
    im_c = im_c .*repmat(imresize(mask,[size(im_c,1) size(im_c,2)]),1,1,3);
elseif ~internet
    % choose some spherical harmonic coefficients for image synthesis
    sh_coeff = [0.0 0.6 0.5 -1.3]/1.2;
    x = sh_coeff(2);   y = sh_coeff(3);   z = -sh_coeff(4);
    A_gt = atan2d(x,z);    E_gt = atan2d(y,z);
    
    Rpose = makehgtform('yrotate',deg2rad(0));
    % Choose image resolution
    imsize = [240 240];
    [im,im_c,z_gt,scales]=read_render_USF(impath,Rpose,imsize);
    im_c_gt = im_c;    alb_gt = im;
    
    [n_gt,N_gnd]=normal_from_depth(z_gt);
    if ~is_albedo
        sh_coeff = sh_coeff/2;
        %set albedo to 1
        im_c = im_c*0+1;
    end
    im_c = render_model_noGL(n_gt,sh_coeff,im_c,0);
    %     im_c = imresize(im2double(imread('.\data\testface.jpg')),1);
else
    im_c = im2double(imread([internet_path internet_imgs{i}]));
end

im = rgb2gray(im_c);
if ~is_albedo
    im = im_c(:,:,1);
end
%% Run face tracker
landmarks = stasm_tracker(im_c,talk);
% landmarks = stasm_tracker(im,talk);
labels = mark_regions(landmarks,im_c);
if isempty(landmarks)
    warning('no landmarks detected by face tracker')
end

%% no texture version
% % im_c = render_model_noGL(n_gt,sh_coeff,im_c,0);
% % im = im_c(:,:,1);
% im = rgb2gray(im_c);

%% Compute pose
% Is the face fronto parallel?
restrictive = 0;
cRes = size(im,2)*0+512;
rRes = size(im,1)*0+512;
[Rpose, Scale] = compute_pose_USF(landmarks, talk, im,restrictive,cRes,rRes);

% Rpose = Rpose/R;
%% generate ref depth map
ply_path = '.\data\ref_model.ply';
dataset_path = '..\datasets\USF 3D Face Data\USF Raw 3D Face Data Set\data_files\test';

[dmap_ref, n_ref, N_ref, alb_ref,eye_mask,scalez] = generate_ref_depthmap_USF(Scale,Rpose,im,im_c,dataset_path,talk);
N_ref(isnan(im))=nan;
%     N_gnd(isnan(N_ref))=NaN;
n_ref((isnan(repmat(im,1,1,3)))) = nan;
dmap_ref(isnan(im))=nan;
%     im(isnan(dmap_ref))=  nan;
if ~is_albedo
    alb_ref = alb_ref*0+1;
end
alb_ref(alb_ref>0.1) = 1;
%% estimate lighting

is_ambient = 1;
non_lin = 0;
l_est_amb_lin = estimate_lighting(n_ref, alb_ref*0+1, im,4,talk,is_ambient,non_lin);
x = l_est_amb_lin(2);   y = l_est_amb_lin(3);   z = -l_est_amb_lin(4);
A_est_amb_lin = atan2d(x,z);    E_est_amb_lin = atan2d(y,z);

Rpose = eye(4);
xrange = [-150 150];
yrange = [-150 150];
c4_amb_lin = render_model_general('./data/sphere.ply', l_est_amb_lin, Rpose, 1000, 1000, xrange, yrange, talk);


figure(f)
subplot(1,2,1)
imshow(im)
% title(sprintf('Ground truth\n A: %.0f, E: %.0f',A_gt,E_gt));

subplot(1,2,2)
imshow(c4_amb_lin)
title(sprintf('Ambient, lin\n A: %.0f, E: %.0f',A_est_amb_lin,E_est_amb_lin));

%% Flow computation
% render reference model for flow computations between this rendering and
% input image
% l_est_amb_lin = l_est_amb_lin;
im_rendered = render_model_noGL(n_ref,l_est_amb_lin,alb_ref,0);
alb_ref2 = alb_ref;
if flow
    % number of pyramid levels
    n_levels = 2;
    % number of iterations of flow on each level
    n_iters = 1;
    % do morphing before optical flow?
    morph = 0;
    dmap_ref_old = dmap_ref;
    [alb_ref2,dmap_ref,eye_mask]=compute_flow(im,im_rendered,landmarks,...
        alb_ref,dmap_ref,eye_mask,n_levels,n_iters,morph,1);
end
dmap_ref(isnan(alb_ref2)) = nan;

% if ~is_albedo
%     l_est_amb_lin = sh_coeff;
% end
dmap_ref(isnan(im))=nan;
dmap_ref(dmap_ref==0) = nan;

n_ref = normal_from_depth(dmap_ref);

is_ambient = 1;
non_lin = 0;
l_est_amb_lin = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient,non_lin);
x = l_est_amb_lin(2);   y = l_est_amb_lin(3);   z = -l_est_amb_lin(4);
A_est_amb_lin = atan2d(x,z);    E_est_amb_lin = atan2d(y,z);


%% Optimization
%     im(isnan(im)) = 0;
if combined
    [depth,alb,is_face] = estimate_depth_alb_nonlin(alb_ref2*255,im*255,dmap_ref,l_est_amb_lin,...
        lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,is_alb_dz,z_gt,algo);
    %         alb = alb*255;
    figure;imshow(alb/255);
    alb_render = alb/255;
    title('estimated albedo')
elseif is_alb_opt
    im_inp = im;
    if raw
        im_inp = im.^(1/2.2);
    end
    [depth,alb,is_face] = estimate_depth_nonlin(alb_ref2*255,im_inp*255,dmap_ref,l_est_amb_lin,...
        lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,is_dz_depth,lambda_dz,[],algo,0);
    figure;imshow(alb/255);
    alb_render = alb/255;
    title('estimated albedo')
else
    depth = estimate_depth_nonlin(alb_ref2*255,im*255,dmap_ref,l_est_amb_lin,...
        lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,z_gt,algo,0);
    alb = im_c*255;
    alb_render = alb_ref2*0+1;
end
is_face = remove_bad_boundary(is_face)>0;

%% Display results
alb_render = max(alb_render,0);
figure;
im_c(im_c<0) = 0;
im_inp_c = im_c;
if raw
    im_inp_c = im_c.^(1/2.2);
end
depth_s=surf(depth,im_inp_c,'edgealpha',0,'facecolor','interp');axis equal
colormap 'gray';
phong.shading(depth_s);
title('Estimated Depth')
if exist('z_gt')
    figure;
    depth_s=surf(z_gt,im_c,'edgealpha',0,'facecolor','interp');axis equal
    colormap 'gray';
    phong.shading(depth_s);
    title('Ground truth')
end
[n_new,N_ref_new] =normal_from_depth( depth );
p_new = n_new(:,:,1).*N_ref_new;
q_new = n_new(:,:,2).*N_ref_new;
% 
% save('state_aunty.mat');
% sc = 0.6;
% im_c = imresize(im_c,sc);
% depth = imresize(depth,sc)*sc;
% labels = imresize(labels,sc,'nearest');
% dmap_ref = imresize(dmap_ref,sc);
% eye_mask = imresize(eye_mask,sc,'nearest');
% is_face = ~isnan(depth);
% is_face = remove_bad_boundary(is_face);
% scalez = scalez/sc;
im_c(im_c<0) = 0;
if raw
    inp_im_c = im_c;
else
    inp_im_c = im_c.^(2.2);
end
[z_o,alb_o,L_sh_o]=intrinsic_decomposition(inp_im_c,depth,labels,(~isnan(dmap_ref).*eye_mask)>0,...
    scalez*1000);

alb_o(alb_o<0)=0;
figure;
depth_s=surf(z_o,alb_o.^(1/2.2)*2,'edgealpha',0,'facecolor','interp');axis equal
phong.shading(depth_s);


im2=imadjust(im_c,[0 1],[0 1],0.7);
chrom = reshape(im2,[],3);
chrom = chrom./repmat(sqrt(sum(chrom.^2,2)),1,3);
chrom = reshape(chrom,size(im2));
figure;

if is_alb_opt
    im_target = im;
    im_target(isnan(N_ref))= 0;
else
    im_target = render_model_noGL(n_gt,l_est_amb_lin,alb_render,0);
    im_target(isnan(N_ref))= 0;
end
subplot(2,2,1);imshow(im_target)
title('Target Image')
subplot(2,2,2);imshow(alb_render)
title('Albedo')
subplot(2,2,3);imshow(render_model_noGL(n_new,l_est_amb_lin,alb_render*0+1,0));
title('Shape')
eye_mask(isnan(eye_mask))=0;
eye_mask = double(eye_mask>0);

subplot(2,2,4);imshow(render_model_noGL(n_new,l_est_amb_lin,alb_render,0));
title('Rendered')
offset = mean(depth(~(isnan(depth) | isnan(z_gt) )))-mean(z_gt(~(isnan(depth) | isnan(z_gt))));
depth2 = depth - offset;
figure;
subplot(1,2,1)
z_error = (depth2-z_gt)*scalez*100;
z_error(isnan(z_error)) = 0;
imagesc(abs(z_error));
title('|z_{est}-z_{ground truth}| _{(cm)}')
if size(alb,3)>1
    alb = rgb2gray(alb/255)*255;
end
im_diff = alb./N_ref_new.*(l_est_amb_lin(2)*p_new+l_est_amb_lin(3)*q_new-l_est_amb_lin(4))-im*255;
im_diff(isnan(im_diff))=0;
im_diff(eye_mask==0) = 0;
subplot(1,2,2);
imagesc(abs(im_diff)*100/255);
title('error in rendered image (Scale of 0-100)')

drawnow
% depth = estimate_depth(N_ref,alb_ref,im,dmap_ref,l,50);

%% Re-lighting
count = count+1;
is_ambient = 1;
non_lin = 0;
l_new = estimate_lighting(n_new,alb_ref2,im,4,0,is_ambient,non_lin);

frameRate = 10;
speed = 5;
span = 0.5;
if raw
    do_relighting_linear( im_c,l_new,n_new,frameRate,speed,span)
else
    do_relighting( im_c,l_new,n_new,frameRate,speed,span)
end