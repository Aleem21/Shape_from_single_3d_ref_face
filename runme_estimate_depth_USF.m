clear variables
close all
import plot3D_helper.label_axis
fusion_path = '.\data\fusion.jpg';

%% set initial variables
folder_path = '.\data\USF_images\';
talk = 0;
impaths = {'03721c15.eko'};
n = numel(impaths);
f=figure;hold on
count = 1;


for i=1:1
    impath = [folder_path impaths{i}];
%% make image
sh_coeff = [0 0.5 0.5 -1.3];
x = sh_coeff(2);   y = sh_coeff(3);   z = -sh_coeff(4);
A_gt = atan2d(x,z);    E_gt = atan2d(y,z);

Rpose = makehgtform('yrotate',deg2rad(0));
[im,im_c,z_gt]=read_render_USF(impath,Rpose);
[n_gt,N_gnd]=normal_from_depth(z_gt);
im_c = render_model_noGL(n_gt,sh_coeff/2,im_c,0);
im = rgb2gray(im_c);
%% Run face tracker
landmarks = stasm_tracker(im,talk);
% landmarks = stasm_tracker(im,talk);

if isempty(landmarks)
    continue;
end

%% no texture version
im_c = render_model_noGL(n_gt,sh_coeff/2,im_c*0+1,0);
im = im_c(:,:,1);
% im = rgb2gray(im_c);

%% Compute pose
restrictive = 0;
[Rpose, Scale] = compute_pose_USF(landmarks, talk, im,restrictive);

% Rpose = Rpose/R;
%% generate ref depth map
ply_path = '.\data\ref_model.ply';
% ply_path ='.\data\sphere.ply';
cRes = size(im,2);
rRes = size(im,1);
[dmap_ref, n_ref, N_ref, alb_ref] = generate_ref_depthmap_USF(Scale,Rpose,im,im_c,talk);
% [dmap_ref, n_ref] = generate_ref_depthmap(ply_path,Scale, talk, 1000, 1000, Rpose,im);


alb_ref = alb_ref*0+1;
% alb_ref = dmap_ref*0+1;

%% estimate lighting
is_ambient = 1;
non_lin = 0;
l_est_amb_lin = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient,non_lin);
x = l_est_amb_lin(2);   y = l_est_amb_lin(3);   z = -l_est_amb_lin(4);
A_est_amb_lin = atan2d(x,z);    E_est_amb_lin = atan2d(y,z);

is_ambient = 0;
l_est = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient);
x = l_est(2);   y = l_est(3);   z = -l_est(4);
A_est = atan2d(x,z);    E_est = atan2d(y,z);


% is_ambient = 1;
% non_lin = 1;
% l_est_amb = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient,non_lin);
% x = l_est_amb(2);   y = l_est_amb(3);   z = -l_est_amb(4);
% A_est_amb = atan2d(x,z);    E_est_amb = atan2d(y,z);


Rpose = eye(4);
xrange = [-150 150];
yrange = [-150 150];
c4 = render_model_general('./data/sphere.ply', l_est, Rpose, 1000, 1000, xrange, yrange, talk);
c4_amb_lin = render_model_general('./data/sphere.ply', l_est_amb_lin, Rpose, 1000, 1000, xrange, yrange, talk);
% c4_amb = render_model_general('./data/sphere.ply', l_est_amb, Rpose, 1000, 1000, xrange, yrange, talk);


figure(f)
subplot(1,2,1)
imshow(im_c)
title(sprintf('Ground truth\n A: %.0f, E: %.0f',A_gt,E_gt));

% subplot(1,2,2)
% imshow(c4_amb)
% title(sprintf('Recovered\n A: %.0f, E: %.0f',A_est_amb_lin,E_est_amb_lin));


% subplot(4,min(n,3),count)
% imshow(im_c)
% title(sprintf('Ground truth\n A: %.0f, E: %.0f',A_gt,E_gt));



% subplot(4,min(n,3),count+min(n,3))
% imshow(c4_amb_lin)
% title(sprintf('Ambient, linear\n A: %.0f, E: %.0f',A_est_amb_lin,E_est_amb_lin));
% 
% 
% subplot(4,min(n,3),count+2*min(n,3))
% imshow(c4)
% title(sprintf('No Ambient, linear\n A: %.0f, E: %.0f',A_est,E_est));
% 
% subplot(4,min(n,3),count+3*min(n,3))
% imshow(c4_amb)
% title(sprintf('Ambient, nonlinear\n A: %.0f, E: %.0f',A_est_amb,E_est_amb));



% im2 = render_model_noGL(n_ref,l_est,alb_ref,talk);
% subplot(3,n,count+2*n)
% imshow(c9)
N_gnd(isnan(N_ref))=NaN;
N_ref = N_gnd;
n_ref = n_gt;
dmap_ref = z_gt-max(z_gt(:));
talk = 1;
l_est = sh_coeff/2;
N_ref_cur = N_ref;
for j = 1:1
    depth = estimate_depth(N_ref_cur,alb_ref,im,dmap_ref,l_est,10,'laplac');
    
    [ ~,N_ref2 ] = normal_from_depth( depth );
        N_ref_cur = (N_ref_cur+N_ref2)/2;
%     N_ref = N_ref2;
    
    depths{j} = depth;
end
figure; surf(depth);axis equal
[n_new,N_ref_new] =normal_from_depth( depth );
p = n_ref(:,:,1).*N_ref;
q = n_ref(:,:,2).*N_ref;
p_new = n_new(:,:,1).*N_ref_new;
q_new = n_new(:,:,2).*N_ref_new;
figure;
target = im;
target(isnan(N_ref))= 0;
subplot(2,2,1);imshow(target)
title('Target Image')
subplot(2,2,2);imshow(alb_ref./N_ref.*(l_est(2)*p+l_est(3)*q-l_est(4)))
title('p_{ref}/N_{ref}')
subplot(2,2,3);imshow(alb_ref./N_ref.*(l_est(2)*p_new+l_est(3)*q_new-l_est(4)))
title('p_{est}/N_{ref}')
subplot(2,2,4);imshow(alb_ref./N_ref_new.*(l_est(2)*p_new+l_est(3)*q_new-l_est(4)))
title('p_{est}/N_{est}')
offset = mean(depth(~(isnan(depth) | isnan(z_gt) )))-mean(z_gt(~(isnan(depth) | isnan(z_gt))));
depth2 = depth - offset;
figure;imagesc(depth2-z_gt);
title('z_{est}-z_{ground truth}')
im_diff = alb_ref./N_ref_new.*(l_est(2)*p_new+l_est(3)*q_new-l_est(4))-im;
im_diff(isnan(im_diff))=0;
figure;imagesc(im_diff);
title('error in rendered image')

drawnow
% depth = estimate_depth(N_ref,alb_ref,im,dmap_ref,l,50);

if mod(count,3)==0 && i<n
    f=figure;
    count = 0;
end
count = count+1;
end