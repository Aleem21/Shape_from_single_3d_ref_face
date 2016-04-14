clear variables
close all
import plot3D_helper.label_axis
fusion_path = '.\data\fusion.jpg';
addpath ('.\SHT')
addpath(genpath('./Basel Face Model'))
%% set initial variables
% optimization weights
lambda1 = 5;
lambda2 = 5;

% number of maximum iterations of depth optimization
max_iter = 500;

% do you want to put in albedo into input image synthesis (keep it 1,
% option was created for early experiments to create the problems easier
% with no albedo complications
is_albedo = 1;

% Keep fixed. Legacy options
lambda_bound = 1;
is_dz_depth = 0;
lambda_dz = 0;
lambda_reg2 = 5;
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
use_model = 1;
% Use some real image as input rather than synthesizing one from dataset?
internet = 1;
internet_path = '.\data\internet\';
internet_imgs = {'person16.jpg'};
fusion_path = '.\data\fusion.jpg';
yale = 0;
raw = 0;
raw_path = '.\data\RAW\';
raw_imgs = {'mahmood2.cr2'};
% % raw_imgs = {'DSC_0190.NEF'};

% USF images path
folder_path = '.\data\USF_images\';

% change this for different USF images
impaths = {'03671c16.eko'};
n = numel(impaths);
f=figure;hold on
count = 1;

% talk = 4 - fully verbose
% talk = 0 - least talking and result displaying.
% talk is not fully consistent in all the code. Needs improvement.
talk = 0;
RGB_to_XYZ = [  0.4124564 0.3575761 0.1804375;
    0.2126729 0.7151522 0.0721750;
    0.0193339 0.1191920 0.9503041];
XYZ_to_RGB = inv(RGB_to_XYZ);
i=1;
impath = [folder_path impaths{i}];
%% make image
if raw && ~ yale
    [im_c,im] = imread_raw([raw_path raw_imgs{i}]);
%     im_c = im_c.^(1/1.5);im = im.^(1/1.5);
    %     im_c = im_c/nanmax(im_c(:));
    %     im = im/nanmax(im(:));
    %     im = rgb2gray(rgb);
    %     im_c = rgb.^(2.2);
    %mask = double(imread([raw_path raw_imgs{i}(1:end-4) '_mask.tiff'])>0);
    %mask = mask(:,:,1);
    sc = 450/max(size(im));
    im_c = imresize(im_c,sc);
    im = imresize(im,sc);
    im_c_full = im_c;
    im_c_full(im_c_full<0) = 0;
    im_full = im;
    %im_c = im_c .*repmat(imresize(mask,[size(im_c,1) size(im_c,2)]),1,1,3);
    %im = im.*imresize(mask,[size(im_c,1) size(im_c,2)]);
    im(im<0) = 0;
    im_c(im_c<0) = 0;
    im_c(repmat(im<1e-4,1,1,3)) = nan;
    im(im<1e-4)=nan;
    lin_xyz = correct(im_c,RGB_to_XYZ);
    %lin_xyz = lin_xyz/nanmax(lin_xyz(:))*0.5;
elseif ~internet && ~yale
    use_model = 1;
    % choose some spherical harmonic coefficients for image synthesis
    %     sh_coeff = [0.3 0.6 0.5 -1.3]/1.2;
    sh_coeff = [1 0 0 0];
    x = sh_coeff(2);   y = sh_coeff(3);   z = -sh_coeff(4);
    A_gt = atan2d(x,z);    E_gt = atan2d(y,z);
    
    Rpose = makehgtform('yrotate',deg2rad(-0))*...
        makehgtform('translate',[0,-0.25,0])*...
        makehgtform('scale',1.0);
    % Choose image resolution
    sc = 1;
    imsize = [480 360]/sc;
    dr = ceil((610-480)/2/sc);
    dc = ceil((458-360)/2/sc);
    is_full = 1;
    [~,im_c_full,z_gt_f]=read_render_USF(impath,Rpose,imsize*1.4/1.1,is_full);
    [n_gt,N_gnd]=normal_from_depth(z_gt_f);
    im_c_full = render_model_noGL(n_gt,[1 0 0 0],im_c_full,0);
    im_c_full_small = im_c_full(dr+[0:floor(480/sc)-1],dc+[0:ceil(360/sc)-1],:);
    is_full = 0;
    [im,im_c,z_gt,scales]=read_render_USF(impath,Rpose,imsize*1.4/1.1,is_full);
    z_gt = z_gt(dr+[0:floor(480/sc)-1],dc+[0:ceil(360/sc)-1]);
    im_c = im_c(dr+[0:floor(480/sc)-1],dc+[0:ceil(360/sc)-1],:);
    im = rgb2gray(im_c);
    im_c_gt = im_c;    alb_gt = rgb2gray(im_c);
    landmarks = stasm_tracker(im_c_full.^(1/2.2),talk);
    landmarks = landmarks - repmat([dc;dr],1,size(landmarks,2));
    labels = mark_regions(landmarks,im_c);
    lin_xyz = rgb2xyz(im_c);
    lin_rgb = correct(lin_xyz,XYZ_to_RGB);
    mask = ones(size(im))>0;
    [D,S,chrom_p] = seperate_specular(lin_rgb,mask,labels);
    lin_xyz = correct(D,RGB_to_XYZ);
    im_c = xyz2rgb(lin_xyz);
%     im_c = D;
    %         D = xyz2rgb(D);
    %         D = rgb2xyz(D);
    %     lin_xyz = lin_xyz/max(max(lin_xyz(:,:,2)))*0.2;
    %     im_c = D;
    % 	im = lin_xyz(:,:,2);
    
    
    
    %     im_c = imresize(im_c,0.8,'nearest');
    %     z_gt = imresize(z_gt,0.8,'nearest');
%     g = fspecial('gaussian',5,1.5);
%     z_gt = conv2(z_gt,g,'same');
    [n_gt,N_gnd]=normal_from_depth(z_gt);
    if ~is_albedo
        sh_coeff = sh_coeff/2;
        %set albedo to 1
        im_c = im_c*0+1;
    end
    %     lin_xyz = rgb2xyz(im_c);
    %     lin_rgb = correct(lin_xyz,XYZ_to_RGB);
    %             lin_rgb = im_c;
%     n_synth = n_gt;
%     n_synth(:,1:end/2,:)= n_gt(:,end:-1:end/2+1,:);
    im_synth = im_c;
    im_synth(:,1:end/2,:)= im_c(:,end:-1:end/2+1,:);
%     im_synth(:,end:-1:end/2+1,:)= im_c_gt(:,1:end/2,:);

	alb_gray = rgb2gray(im_synth);
%     alb_gray(:) = 1;
    im_c = zeros(size(n_gt(:,:,1)));
    for i=1:3
        sh_coeff = [0;rand(2,1)*2-1;-2];
        sh_coeff = sh_coeff/norm(sh_coeff(2:4))*0.33;
        im_c = im_c + max(render_model_noGL(n_gt,sh_coeff,alb_gray/max(alb_gray(:)),0),0);
    end
    %      sh_coeff = [1];
    %     lin_rgb(lin_rgb>0) = 1;
    
    im_c(im_c==1)=nan;
    lin_xyz = correct(im_c,RGB_to_XYZ);
    lin_xyz(lin_xyz==1)=nan;
    lin_xyz = lin_xyz;
    %     scl = nanmax(lin_xyz(:))/0.5;
    % 	lin_xyz = lin_xyz/scl;
    %     sh_coeff = sh_coeff/scl;
%     im = lin_xyz(:,:,2);
    im = im_c;
    im = real(im);
    
    im_c_full = im_c_full.^(2.2);
    %     im_c = imresize(im2double(imread('.\data\testface.jpg')),1);
elseif ~yale
    
    im_c = im2double(imread([internet_path internet_imgs{i}])).^2.2;
%     gamma = blindgamma(1,0.4:0.1:2.4,im_c);
%     im_c = imresize(im_c,0.6);
    
    %gauss = fspecial('gaussian',3,1.5);
    %im_c = imfilter(im_c,gauss);
    im_c(im_c<0) =0;
    landmarks = stasm_tracker(im_c,talk);
    labels = mark_regions(landmarks,im_c);
    mask = ones(size(im_c(:,:,1)))>0;
    [D,S]=seperate_specular(correct(rgb2xyz(im_c),XYZ_to_RGB),mask,labels);
    D = rgb2xyz(im_c);
    %     lin_xyz = rgb2xyz(im_c);
    lin_xyz= correct(D,RGB_to_XYZ);
    %lin_xyz = imfilter(lin_xyz,gauss);
    %     im_c(repmat(lin_xyz(:,:,2)==0,1,1,3)) = nan;
    %     lin_xyz = rgb2xyz(im_c);
    %     im_c = im_c.^(2.4);
    %     im = rgb2gray(im_c);
    %     lin_xyz = lin_xyz;
    im = lin_xyz(:,:,2);
%     im_c = im_c.^(2.2);
    im = rgb2gray(im_c); %% Watch out for this non linear thing!
%         im_c_full = im_c.^(2.2);
    % im_c_full = imadjust(im_c,[0 1],[0 1],2);
    im_c_full = im_c;
else
    folder_path = '.\data\yaleB_l_test\';
    impaths = {'yaleB08_P00A+000E+00.pgm','yaleB04_P00A+000E+20.pgm',...
        'yaleB08_P00A+000E+20.pgm','yaleB01_P00A+000E+90.pgm'...
        ,'yaleB01_P00A+025E+00.pgm','yaleB01_P00A+050E+00.pgm',...
        'yaleB01_P00A+070E+45.pgm','yaleB27_P00A+000E+00.pgm',...
        'yaleB29_P00A+000E+00.pgm','yaleB29_P00A+000E+20.pgm'};
    impath = [folder_path impaths{end}];
    
    [im,im_c] = read_image(impath);
    im_c = repmat(im_c,1,1,3);
    fused = fuse_image_YaleB(fusion_path, im);
    im_c_full = im_c;
    %im = im/max(im(:))*0.5;
    use_model = 0;
end

if ~is_albedo
    im = im_c(:,:,1);
end
%% Run face tracker
if ~raw && ~internet && ~yale
    landmarks = stasm_tracker(im_c_full.^(1/2.2),talk);
    landmarks = landmarks - repmat([dc;dr],1,size(landmarks,2));
elseif yale
    landmarks = stasm_tracker_YaleB(fused,talk);
%     im = imresize(im,2);
%     im_c = imresize(im_c,2);
%     im_c_full = im_c;
else
    landmarks = stasm_tracker(imresize(im_c,2),talk);
    landmarks = floor(landmarks/2);
end
% landmarks = stasm_tracker(im,talk);
if isempty(landmarks)
    warning('no landmarks detected by face tracker')
end
labels = mark_regions(landmarks,im_c);
%% spec removal
%     mask = ones(size(im))>0;
%     [D,S,chrom_p] = seperate_specular(correct(lin_xyz,XYZ_to_RGB),mask,labels);
%     lin_xyz = correct(D,RGB_to_XYZ);
% %         D = xyz2rgb(D);
% %         D = rgb2xyz(D);
% %     lin_xyz = lin_xyz/max(max(lin_xyz(:,:,2)))*0.2;
%     im_c = D;
% 	im = lin_xyz(:,:,2);


%% synth new model
%     load('./Basel Face Model/01_MorphableModel.mat');
%     load('./Basel Face Model/02_scans_matlab/00002_20061015_00448_neutral_face05.mat')
%     FV.Vertices = reshape(shapeMU+5e5*shapePC(:,11),3,[])';
%     FV.FaceVertexCData = abs(reshape(texMU+5e3*texPC(:,9),3,[])'/255);
%
%     r = makehgtform( 'yrotate',-0.0);
%     r = r(1:3,1:3);
%     FV.Vertices = FV.Vertices*r';
%     FV.Faces = tl;
%     FV.FaceColor = 'interp';
% %     FV.FaceVertexCData = tex'/255;
%
%     xrange = [1 240];
%     yrange = [1 240];
%     pts= double(FV.Vertices')/9e2 + 120;
%     pts(1,:) = -(pts(1,:) - mean(xrange)) + mean(xrange);
%     pts(3,:) = -(pts(3,:) - max(pts(3,:)));
%     tri = FV.Faces'-1;
%     clr = double(FV.FaceVertexCData');
%     [depth_map,albedo] = computer_depth_USF( pts,tri,clr,xrange,yrange,[],0 );
%     figure;imshow(albedo)
% %         im = rgb2xyz(im_c);
% %     im = im(:,:,2);
%     albedo_xyz = rgb2xyz(albedo);
%     albedo_c = correct(albedo_xyz,XYZ_to_RGB);
%     albedo_g = albedo_xyz(:,:,2);
%     albedo_g = rgb2gray(albedo);
% %     depth_map(177,103) = 23.56;
%     n_gt = normal_from_depth(depth_map);
%     l_synth = [0.1;0.05;-0.02;-0.1]*2;
%     im_c = render_model_noGL(n_gt,l_synth,albedo,1);
%     albedo_g = rgb2gray(albedo);
%     im_c_full = im_c;
%     im_xyz = correct(im_c,RGB_to_XYZ);
%     im = im_xyz(:,:,2);
%     im = rgb2gray(im_c);
%     %     im_xyz = correct(im_c,RGB_to_XYZ);
%     %     im = im_xyz(:,:,2)*255;
%     %     im_c = im_c*255;
%     %     im_c = correct(im_xyz,XYZ_to_RGB);
%     load('mask')
%         clr = double(repmat(mask',3,1)==0);
%     [~,mask] = computer_depth_USF( pts,tri,clr,xrange,yrange,[],0 );
%     im(mask(:,:,1)==0) = nan;
%     im_c(mask==0) = nan;
%     landmarks = stasm_tracker(im_c_full,talk);
%     labels = mark_regions(landmarks,im_c_full);


%% no texture version
% % im_c = render_model_noGL(n_gt,sh_coeff,im_c,0);
% % im = im_c(:,:,1);
% im = rgb2gray(im_c);

%% Compute pose
% Is the face fronto parallel?
restrictive = 1;
cRes = size(im,2)*0+512;
rRes = size(im,1)*0+512;
[Rpose, Scale] = compute_pose_USF(landmarks, talk, im,restrictive,cRes,rRes);

% Rpose = Rpose/R;
%% generate ref depth map
ply_path = '.\data\ref_model.ply';
dataset_path = '..\datasets\USF 3D Face Data\USF Raw 3D Face Data Set\data_files\test';
warning off
if use_model
    [dmap_ref, n_ref, N_ref, alb_ref,eye_mask,scalez,~,model] = generate_ref_depthmap_USF(Scale,Rpose,rgb2gray(im_c_full),im_c_full.^(1/2.2),dataset_path,use_model,talk);
    if ~internet && ~yale && ~raw
        dmap_ref = dmap_ref(dr+[0:floor(480/sc)-1],dc+[0:ceil(360/sc)-1]);
        n_ref = n_ref(dr+[0:floor(480/sc)-1],dc+[0:ceil(360/sc)-1],:);
        alb_ref = alb_ref(dr+[0:floor(480/sc)-1],dc+[0:ceil(360/sc)-1],:);
        eye_mask = eye_mask(dr+[0:floor(480/sc)-1],dc+[0:ceil(360/sc)-1]);
    end
else
    [dmap_ref, n_ref, N_ref, alb_ref,eye_mask,scalez] = generate_ref_depthmap_USF(Scale,Rpose,im,im_c,dataset_path,use_model,talk);
end
im_calib = im;
im_calib(isnan(dmap_ref)) = nan;
% im = im/nanmedian(im_calib(:))*0.338;
% im_c = im_c/nanmedian(im_calib(:))*0.338;
% dmap_ref(154,120) = dmap_ref(154,121);
warning on
labels2 = mark_regions(landmarks,im_c,1);
eye_mask = eye_mask==1;
% eye_mask(:) = 1;
% dmap_ref = z_gt;
% n_ref = n_gt;
% N_ref = N_gnd;
eye_mask = double(eye_mask);
im(isnan(dmap_ref)) = nan;
dmap_ref(isnan(im))=nan;
alb_ref(alb_ref==0)=nan;
% alb_ref(:) = 1;
N_ref(isnan(alb_ref))=nan;
N_ref(isnan(im))=nan;
%     N_gnd(isnan(N_ref))=NaN;
n_ref((isnan(repmat(im,1,1,3)))) = nan;
dmap_ref(isnan(im))=nan;
n_ref((isnan(repmat(dmap_ref,1,1,3)))) = nan;
%     im(isnan(dmap_ref))=  nan;
if ~is_albedo
    alb_ref = alb_ref*0+1;
end
% alb_ref =alb_ref*0+1;
% alb_ref(alb_ref>0.1) = 1;
%% synthetic magic
% l_synth = [0.05;0.01;-0.02;-0.1]*5;
% dmap_ref = z_gt;
% % dmap_ref = imresize(dmap_ref,0.5,'nearest');
% n_ref = normal_from_depth(dmap_ref);
% % alb_ref = imresize(alb_ref,0.5,'nearest');
% n_synth = normal_from_depth(dmap_ref);
%

% alb_ref(:) = 1;
% im = render_model_noGL(n_synth,l_synth,alb_ref,0);
%
% % labels = imresize(labels,0.5,'nearest');
% % im_c = imresize(im_c ,0.5,'nearest');

%% estimate lighting
% dmap_ref = dmap_ref*2;
% [n_ref,N_ref] = normal_from_depth(dmap_ref);
is_ambient = 1;
non_lin = 0;
l_est_amb_lin = estimate_lighting(n_ref, alb_ref*255, im*255,4,talk,is_ambient,non_lin,eye_mask);
x = l_est_amb_lin(2);   y = l_est_amb_lin(3);   z = -l_est_amb_lin(4);
A_est_amb_lin = atan2d(x,z);    E_est_amb_lin = atan2d(y,z);
% l_est_amb_lin  = sh_coeff;
im_rendered = render_model_noGL(n_ref,l_est_amb_lin,alb_ref,0);

is_ambient = 1;
l_est_amb_lin2 = estimate_lighting(n_ref, alb_ref*255, im*255,9,talk,is_ambient,non_lin,eye_mask);
im_rendered2 = render_model_noGL(n_ref,l_est_amb_lin2,alb_ref,0);
figure;subplot(1,3,1);imshow(im);subplot(1,3,2);imshow(im_rendered);subplot(1,3,3);imshow(im_rendered2);
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
sh2sph(l_est_amb_lin,[300,300],1);

% is_ambient = 1;
% non_lin = 0;
% l_est_amb_lin2 = estimate_lighting(n_ref, alb_ref2*255, im*255,9,talk,is_ambient,non_lin,eye_mask);
% shadow_ratio = get_shadow(model,dmap_ref,l_shadow);
%[s_ratio,im_rendered_s,im_rendered_ns] = get_shadow_z(dmap_ref,im,l_est_amb_lin);

%alb_ref2 = alb_ref.*s_ratio;
if flow
    % number of pyramid levels
    n_levels = 3;
    % number of iterations of flow on each level
    n_iters = 2;
    % do morphing before optical flow?
    morph = 1;
%         morph = 0;

    dmap_ref_old = dmap_ref;
    im = double(real(im));
    im(im<0) = 0;
    im_rendered = double(real(im_rendered));
    im_rendered(im_rendered<0) = 0;
    if yale
        fused = fuse_image_YaleB(fusion_path, im_rendered);
        yale_landmarks = stasm_tracker_YaleB(fused,0);
        [alb_ref2,dmap_ref2,eye_mask]=compute_flow(im,im_rendered,landmarks,...
            alb_ref,dmap_ref,eye_mask,n_levels,n_iters,morph,1,yale_landmarks);
    elseif raw
        [alb_ref2,dmap_ref2,eye_mask]=compute_flow(im,im_rendered,landmarks,...
            alb_ref,dmap_ref,eye_mask,n_levels,n_iters,morph,1);
    else
        [alb_ref2,dmap_ref2,eye_mask]=compute_flow(im,im_rendered,landmarks,...
            alb_ref,dmap_ref,eye_mask,n_levels,n_iters,morph,1);
    end
else
    dmap_ref2 = dmap_ref;
end



dmap_ref2(isnan(alb_ref2)) = nan;
eye_mask(isnan(eye_mask))= 1;
% eye_mask = eye_mask & labels2~=3 & labels2~=4;
% if ~is_albedo
%     l_est_amb_lin = sh_coeff;
% end
% dmap_ref2 = dmap_ref2 * 1.5;
dmap_ref2(isnan(im))=nan;
dmap_ref2(isnan(alb_ref2)) = nan;
dmap_ref2(dmap_ref2==0) = nan;
n_ref_old = n_ref;

n_ref = normal_from_depth(dmap_ref2);
% n_ref = normal_from_depth(z_gt);
n_ref((isnan(repmat(im,1,1,3)))) = nan;
% alb_ref2 = alb_ref2*0+1;

%% Estimate shadow
%im_rendered = render_model_noGL(n_ref,l_est_amb_lin,alb_ref2,0);
is_ambient = 1;
non_lin = 0;
l_shadow = estimate_lighting(n_ref, alb_ref2*255, im*255,4,talk,is_ambient,non_lin,eye_mask);
im_rendered_s = render_model_noGL(n_ref,l_shadow,alb_ref2,0);
l_shadow9 = estimate_lighting(n_ref,alb_ref2*255, im*255,9,talk,is_ambient,non_lin,eye_mask);
im_rendered_9s = render_model_noGL(n_ref,l_shadow9,alb_ref2,0);

% shadow_ratio = get_shadow(model,dmap_ref,l_shadow);
[s_ratio,shad,nshad] = get_shadow_z(dmap_ref2,im,l_shadow);
l_shadow9 = estimate_lighting(n_ref,alb_ref2.*255.*s_ratio, im*255,9,talk,is_ambient,non_lin,eye_mask);
[s_ratio9, shad9, nshad9] = get_shadow_z(dmap_ref2,im,l_shadow9,0);
l_shadow2 = estimate_lighting(n_ref, alb_ref2.*255.*s_ratio,im*255,9,talk,is_ambient,non_lin,eye_mask);
%%

alb_ref2_9 = alb_ref2.*s_ratio9;
alb_ref2_4 = alb_ref2.*s_ratio;

% alb_ref2  = double(albedo_g);
% im = im.^(1/2.2);
% im = im/10;
% im = im/max(im(:))*0.5;

dmap_ref2(isnan(im))=nan;
gauss = fspecial('gaussian',20,5);
eye_mask = conv2(double(eye_mask),gauss,'same')./conv2(double(eye_mask)*0+1,gauss,'same');

l_est_amb_lin_o = l_est_amb_lin;
is_ambient = 1;
non_lin = 0;
[l_est_amb_lin] = estimate_lighting(n_ref, alb_ref2_9*255, im*255,4,talk,is_ambient,non_lin,eye_mask);
[l_est_amb_lin_ns] = estimate_lighting(n_ref, alb_ref2*255, im*255,4,talk,is_ambient,non_lin,eye_mask);
x = l_est_amb_lin(2);   y = l_est_amb_lin(3);   z = -l_est_amb_lin(4);
A_est_amb_lin = atan2d(x,z);    E_est_amb_lin = atan2d(y,z);

% l_est_amb_lin = l_synth;
%% Optimization
%     im(isnan(im)) = 0;
if combined
    [depth,alb,is_face] = estimate_depth_alb_nonlin(alb_ref2*255,im*255,dmap_ref2*0+255,l_est_amb_lin,...
        lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,is_alb_dz,z_gt,algo);
    %         alb = alb*255;
    figure;imshow(alb/255);
    alb_render = alb/255;
    title('estimated albedo')
elseif is_alb_opt
    im_inp = im;
    %     if raw
    %         im_inp = im.^(1/2.2);
    %     end
    RGB_to_XYZ = [  0.4124564 0.3575761 0.1804375;
        0.2126729 0.7151522 0.0721750;
        0.0193339 0.1191920 0.9503041];
    if raw && ~yale
        %         [D,S,chrom_p] = seperate_specular(im_c,eye_mask>0,labels);
        %         lin_XYZ = correct(D,RGB_to_XYZ);
        %         D = lin_XYZ(:,:,2);
        lin_xyz = correct(im_c,RGB_to_XYZ);
        %         [D,S,chrom_p] = seperate_specular(xyz2rgb(lin_xyz),eye_mask>0,labels);
        [D,S,chrom_p] = seperate_specular(correct(lin_xyz,XYZ_to_RGB),eye_mask>0,labels);
        %         D = rgb2xyz(correct(D);
        D = correct(D,RGB_to_XYZ);
        D = D(:,:,2);
        im = D;
        im = im/max(im(:))*0.5;
        im = im/nanmean(im(:))*0.15;
        n_ref(isnan(repmat(im,1,1,3))) = nan;
        is_ambient = 2;
        l_est_amb_lin = estimate_lighting(n_ref, alb_ref2_9*255, im*255,4,talk,is_ambient,non_lin,eye_mask);
        l_est_amb_lin_ns = estimate_lighting(n_ref, alb_ref2*255, im*255,4,talk,is_ambient,non_lin,eye_mask);
        
        
    else
        
    end
    
    
    lambda_reg2 = 0;
    is_dz_depth = 1;
    lambda_dz = 1;
    light_scale = sqrt(nanmedian(im(:))/0.338);
    reg_scale = 1/sqrt(scalez/3.25e-4);
    lambda1 = 30*light_scale*reg_scale;
    lambda1 = 100;
    % eye_mask(:) = 1;
    lambda_bound = 1;
    max_iter = 20;
    is_l_sh = 0;
    % l_est_amb_lin = l_synth;
    rng(1);
    if ~internet && ~yale && ~ raw
        dmap_ref2(im==0) = nan;
    end
%     [depth] = estimate_depth_nonlin(alb_ref2_9*255,im*255,...
%         dmap_ref2,l_est_amb_lin,...
%         lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask.*s_ratio9,...
%         is_dz_depth,is_l_sh,lambda_dz,lambda_reg2,[],algo,1);
    [depth] = estimate_depth_nonlin(alb_ref2*255,im*255,...
        dmap_ref2,l_est_amb_lin_ns,...
        lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,...
        is_dz_depth,is_l_sh,lambda_dz,lambda_reg2,[],algo,1);
    
    LoG = diag([0 0 1 0 0])-fspecial('gaussian',5,6);
    dz = depth-dmap_ref2;
    highfreq = conv2(dz,LoG,'same');
    
    lambda_dz = 2;
    is_dz_depth = 1;
    lambda1 = 400;
    [depth2] = estimate_depth_nonlin(alb_ref2_9*255,im*255,...
        dmap_ref2,l_est_amb_lin,...
        lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask.*s_ratio9,...
        is_dz_depth,is_l_sh,lambda_dz,lambda_reg2,[],algo,1);
    
%     [depth2] = estimate_depth_nonlin(alb_ref2*255,im*255,...
%         dmap_ref2,l_est_amb_lin_ns,...
%         lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,...
%         is_dz_depth,is_l_sh,lambda_dz,lambda_reg2,[],algo,1);
    
    
    is_face = ~isnan(depth2);
    %     offset = mean(depth(is_face)) - mean(dmap_ref2(is_face));
    %     offset = mean(depth(is_face)) - mean(z_gt(is_face));
    %     depth = depth-offset;
    %     figure;surf(depth);axis equal
    %     figure;imshow(alb/255);
    %     alb_render = alb/255;
    %     title('estimated albedo')
else
    depth = estimate_depth_nonlin(alb_ref2*255,im*255,dmap_ref2,l_est_amb_lin,...
        lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,z_gt,algo,0);
    alb = im_c*255;
    alb_render = alb_ref2*0+1;
end
is_face = remove_bad_boundary(is_face)>0;

%% Display results
faceclr = zeros(size(im_c));
faceclr(:,:,1) = 0.951; faceclr(:,:,2) = 0.6686; faceclr(:,:,3) = 0.3922;
figure;
im_c(im_c<0) = 0;
im_inp_c = im_c;
if ~yale && (raw || internet)
    im_inp_c = xyz2rgb(lin_xyz);
    im_inp_c = xyz2rgb(correct(im_c,RGB_to_XYZ));
    im_inp_c = im_c.^(1/2.2);
else
    im_inp_c = im_c;
end
depth_s=surf(depth,faceclr,'edgealpha',0,'facecolor','interp');axis equal
colormap 'gray';
phong.shading(depth_s);
title('All high')



figure;
im_c(im_c<0) = 0;
im_inp_c = im_c;
if ~yale && (raw || internet)
    im_inp_c = xyz2rgb(lin_xyz);
    im_inp_c = xyz2rgb(correct(im_c,RGB_to_XYZ));
    im_inp_c = im_c.^(1/2.2);
else
    im_inp_c = im_c;
end
depth_s=surf(depth2,faceclr,'edgealpha',0,'facecolor','interp');axis equal
colormap 'gray';
phong.shading(depth_s);
title('All low')


figure;
im_c(im_c<0) = 0;
im_inp_c = im_c;
if ~yale && (raw || internet)
    im_inp_c = xyz2rgb(lin_xyz);
    im_inp_c = xyz2rgb(correct(im_c,RGB_to_XYZ));
    im_inp_c = im_c.^(1/2.2);
else
    im_inp_c = im_c;
end
depth_s=surf(depth2+highfreq*2,faceclr,'edgealpha',0,'facecolor','interp');axis equal
colormap 'gray';
phong.shading(depth_s);
title('High + Low')



figure;
im_c(im_c<0) = 0;
im_inp_c = im_c;
if ~yale && (raw || internet)
    im_inp_c = xyz2rgb(lin_xyz);
    im_inp_c = xyz2rgb(correct(im_c,RGB_to_XYZ));
    im_inp_c = im_c.^(1/2.2);
else
    im_inp_c = im_c;
end
depth_s=surf(depth2+highfreq,im_inp_c,'edgealpha',0,'facecolor','interp');axis equal
colormap 'gray';
%phong.shading(depth_s);
title('High + Low')


if exist('z_gt')
    zgt = z_gt;
    zgt(~is_face) = nan;
    figure;
    depth_s=surf(zgt,im_inp_c*0+0.5,'edgealpha',0,'facecolor','interp');axis equal
    colormap 'gray';
    phong.shading(depth_s);
    title('Ground truth')
end



[n_new,N_ref_new] =normal_from_depth( depth2);
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

% Synthetic image making
[~,S_r,D_r,alpha,P_dr,Yi,Ft_o,E_lm,d_S_n,l_i,Y_N]=render_DS(N_v',w_o,n,rho_s,alb,r,L_sh,is_face,delta_a,delta_s1,chrom);

im_rendered = render_model_noGL(n_new,l_est_amb_lin,alb_ref2,0);



[z_o,alb_o,L_sh_o]=intrinsic_decomposition(inp_im_c,depth,labels,(~isnan(dmap_ref2).*eye_mask)>0,...
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
span = 1;
if raw
    do_relighting_linear( im_c,l_new,n_new,frameRate,speed,span)
else
    do_relighting( im_c,l_est_amb_lin,n_new,frameRate,speed,span)
end

