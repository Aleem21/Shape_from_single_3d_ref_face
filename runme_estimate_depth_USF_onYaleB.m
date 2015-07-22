clear variables
close all
import plot3D_helper.label_axis
fusion_path = '.\data\fusion.jpg';

%% set initial variables
folder_path = '.\data\yaleB_l_test\';
talk = 4;
impaths = {'yaleB08_P00A+000E+00.pgm','yaleB01_P00A+000E+00.pgm', 'yaleB01_P00A+000E+45.pgm','yaleB01_P00A+000E+90.pgm'...
    ,'yaleB01_P00A+025E+00.pgm','yaleB01_P00A+050E+00.pgm','yaleB01_P00A+070E+45.pgm'};
n = numel(impaths);
f=figure;hold on
count = 1;


for i=2:n
    impath = [folder_path impaths{i}];
%% make image
[im,im_c] = read_image(impath);
fused = fuse_image_YaleB(fusion_path, im);

%% Run face tracker
landmarks = stasm_tracker_YaleB(fused,talk);
% landmarks = stasm_tracker(im,talk);

if isempty(landmarks)
    continue;
end

%% Compute pose
restrictive = 1;
[Rpose, Scale] = compute_pose_USF(landmarks, talk, im,restrictive);

% Rpose = Rpose/R;
%% generate ref depth map
ply_path = '.\data\ref_model.ply';
% ply_path ='.\data\sphere.ply';
cRes = size(im,2);
rRes = size(im,1);
[dmap_ref, n_ref, N_ref, alb_ref] = generate_ref_depthmap_USF(Scale,Rpose,im,im_c,talk);
% [dmap_ref, n_ref] = generate_ref_depthmap(ply_path,Scale, talk, 1000, 1000, Rpose,im);



% alb_ref = dmap_ref*0+1;

%% estimate lighting
is_ambient = 1;
non_lin = 0;
l_est_amb_lin = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient,non_lin);
x = l_est_amb_lin(2);   y = l_est_amb_lin(3);   z = -l_est_amb_lin(4);
A_est_amb_lin = atan2d(x,z);    E_est_amb_lin = atan2d(y,z);

is_ambient = 0;
non_lin = 0;
l_est = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient,non_lin);
x = l_est(2);   y = l_est(3);   z = -l_est(4);
A_est = atan2d(x,z);    E_est = atan2d(y,z);

is_ambient = 1;
non_lin = 1;
l_est_amb = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient,non_lin);
x = l_est_amb(2);   y = l_est_amb(3);   z = -l_est_amb(4);
A_est_amb = atan2d(x,z);    E_est_amb = atan2d(y,z);


Rpose = eye(4);
xrange = [-150 150];
yrange = [-150 150];
c4 = render_model_general('./data/sphere.ply', l_est/2, Rpose, 1000, 1000, xrange, yrange, talk);
c4_amb_lin = render_model_general('./data/sphere.ply', l_est_amb_lin/2, Rpose, 1000, 1000, xrange, yrange, talk);
c4_amb = render_model_general('./data/sphere.ply', l_est_amb/2, Rpose, 1000, 1000, xrange, yrange, talk);


figure(f)
subplot(4,min(n,3),count)
imshow(im_c)
title(impath(end-12:end-4))



subplot(4,min(n,3),count+min(n,3))
imshow(c4_amb_lin)
title(sprintf('Ambient, linear\n A: %.0f, E: %.0f',A_est_amb_lin,E_est_amb_lin));


subplot(4,min(n,3),count+2*min(n,3))
imshow(c4)
title(sprintf('No Ambient, linear\n A: %.0f, E: %.0f',A_est,E_est));

subplot(4,min(n,3),count+3*min(n,3))
imshow(c4_amb)
title(sprintf('Ambient, nonlinear\n A: %.0f, E: %.0f',A_est_amb,E_est_amb));



% im2 = render_model_noGL(n_ref,l_est,alb_ref,talk);
% subplot(3,n,count+2*n)
% imshow(c9)
depth = estimate_depth(N_ref,alb_ref,im,dmap_ref,l_est,0.1,'dz');
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