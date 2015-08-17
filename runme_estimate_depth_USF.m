clear variables
close all
import plot3D_helper.label_axis
fusion_path = '.\data\fusion.jpg';

%% set initial variables
lambda1 = 70;
lambda2 = 0.5;
lambda_bound = 1;
max_iter = 50;
is_albedo = 1;
is_alb_opt = 1;
jack = 'on';
combined = 1;
boundary_type = 2;
flow = 1;
folder_path = '.\data\USF_images\';
talk = 0;
impaths = {'03653c16.eko'};
n = numel(impaths);
f=figure;hold on
count = 1;

for i=1:1
    impath = [folder_path impaths{i}];
    %% make image
    sh_coeff = [0 0.6 0.5 -1.3];
    x = sh_coeff(2);   y = sh_coeff(3);   z = -sh_coeff(4);
    A_gt = atan2d(x,z);    E_gt = atan2d(y,z);
    
    Rpose = makehgtform('yrotate',deg2rad(0));
    [im,im_c,z_gt,scales]=read_render_USF(impath,Rpose,[240 240]);
    alb_gt = im;
    %     z_gt(:,[1:40 201:end]) = [];
    %     z_gt(181:end,:) = [];
    %     im_c(:,[1:40 201:end],:) = [];
    %     im_c(:,201:end,:) = [];
    %     im_c(181:end,:,:) = [];
    [n_gt,N_gnd]=normal_from_depth(z_gt);
    if ~is_albedo
        sh_coeff = sh_coeff/2;
        %set albedo to 1
        im_c = im_c*0+1;
    end
    im_c = render_model_noGL(n_gt,sh_coeff,im_c,0);
    %     im_c = imresize(im2double(imread('.\data\l_test\scarlet.jpg')),0.3);
    im = rgb2gray(im_c);
    if ~is_albedo
        im = im_c(:,:,1);
    end
    %% Run face tracker
    landmarks = stasm_tracker(im,talk);
    % landmarks = stasm_tracker(im,talk);
    
    if isempty(landmarks)
        continue;
    end
    
    %% no texture version
    % % im_c = render_model_noGL(n_gt,sh_coeff,im_c,0);
    % % im = im_c(:,:,1);
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
    [dmap_ref, n_ref, N_ref, alb_ref,eye_mask,scalez] = generate_ref_depthmap_USF(Scale,Rpose,im,im_c,talk);
    % [dmap_ref, n_ref] = generate_ref_depthmap(ply_path,Scale, talk, 1000, 1000, Rpose,im);
    
    N_ref(isnan(im))=nan;
    %     N_gnd(isnan(N_ref))=NaN;
    n_ref((isnan(repmat(im,1,1,3)))) = nan;
    dmap_ref(isnan(im))=nan;
    im(isnan(dmap_ref))=  nan;
    if ~is_albedo
        alb_ref = alb_ref*0+1;
    end
    
    %% estimate lighting
    
    is_ambient = 1;
    non_lin = 0;
    l_est_amb_lin = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient,non_lin);
    x = l_est_amb_lin(2);   y = l_est_amb_lin(3);   z = -l_est_amb_lin(4);
    A_est_amb_lin = atan2d(x,z);    E_est_amb_lin = atan2d(y,z);
    
    is_ambient = 0;
    non_lin = 0;
    l_est_nonamb_lin = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient);
    x = l_est_nonamb_lin(2);   y = l_est_nonamb_lin(3);   z = -l_est_nonamb_lin(4);
    A_est_nonamb_lin = atan2d(x,z);    E_est_nonamb_lin = atan2d(y,z);
    
    
    is_ambient = 1;
    non_lin = 1;
    l_est_amb_nonlin = estimate_lighting(n_ref, alb_ref, im,4,talk,is_ambient,non_lin);
    x = l_est_amb_nonlin(2);   y = l_est_amb_nonlin(3);   z = -l_est_amb_nonlin(4);
    A_est_amb_nonlin = atan2d(x,z);    E_est_amb_nonlin = atan2d(y,z);
    
    
    Rpose = eye(4);
    xrange = [-150 150];
    yrange = [-150 150];
    c4_nonamb_lin = render_model_general('./data/sphere.ply', l_est_nonamb_lin, Rpose, 1000, 1000, xrange, yrange, talk);
    c4_amb_lin = render_model_general('./data/sphere.ply', l_est_amb_lin, Rpose, 1000, 1000, xrange, yrange, talk);
    c4_amb_nonlin = render_model_general('./data/sphere.ply', l_est_amb_nonlin, Rpose, 1000, 1000, xrange, yrange, talk);
    
    
    figure(f)
    subplot(2,2,1)
    imshow(im)
    title(sprintf('Ground truth\n A: %.0f, E: %.0f',A_gt,E_gt));
    
    subplot(2,2,2)
    imshow(c4_amb_lin/2)
    title(sprintf('Ambient, lin\n A: %.0f, E: %.0f',A_est_amb_lin,E_est_amb_lin));
    
    
    subplot(2,2,3)
    imshow(c4_nonamb_lin/2)
    title(sprintf('Non-ambient, lin\n A: %.0f, E: %.0f',A_est_nonamb_lin,E_est_nonamb_lin));
    
    
    
    subplot(2,2,4)
    imshow(c4_amb_nonlin/2)
    title(sprintf('Ambient, non-linear\n A: %.0f, E: %.0f',A_est_amb_nonlin,E_est_amb_nonlin));
    
    
    im_rendered = render_model_noGL(n_ref,l_est_nonamb_lin,alb_ref,0);
    if flow
        n_levels = 10;
        n_iters = 1;
        morph = 1;
        [alb_ref2,dmap_ref,eye_mask]=compute_flow(im,im_rendered,...
            landmarks,alb_ref,dmap_ref,eye_mask,n_levels,n_iters,morph,1);
    end
    talk = 1;
    if ~is_albedo
        l_est_nonamb_lin = sh_coeff;
    end
    dmap_ref(isnan(im))=nan;
    %     im(isnan(im)) = 0;
    if combined
        [depth,alb] = estimate_depth_alb_nonlin(alb_ref2,im,dmap_ref,l_est_nonamb_lin,...
            lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,z_gt);
        alb = alb*255;
        figure;imshow(alb/255);
        alb_render = alb/255;
        title('estimated albedo')
    elseif is_alb_opt
        [depth,alb] = estimate_depth_nonlin(alb_ref2*255,im*255,dmap_ref,l_est_nonamb_lin,...
            lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,z_gt);
        figure;imshow(alb/255);
        alb_render = alb/255;
        title('estimated albedo')
    else
        depth = estimate_depth_nonlin(alb_ref2*255,im*255,dmap_ref,l_est_nonamb_lin,...
            lambda1,lambda2,lambda_bound,max_iter,boundary_type,jack,eye_mask,z_gt);
        alb = im_c*255;
        alb_render = alb_ref2*0+1;
    end
    
    figure;
    depth_s=surf(depth,im_c,'edgealpha',0,'facecolor','interp');axis equal
    colormap 'gray';
    phong.shading(depth_s);
    title('Estimated Depth')
    figure;
    depth_s=surf(z_gt,im_c,'edgealpha',0,'facecolor','interp');axis equal
    colormap 'gray';
    phong.shading(depth_s);
    title('Ground truth')
    [n_new,N_ref_new] =normal_from_depth( depth );
    p_new = n_new(:,:,1).*N_ref_new;
    q_new = n_new(:,:,2).*N_ref_new;
    figure;
    im_target = render_model_noGL(n_gt,l_est_nonamb_lin,alb_render,0);
    im_target(isnan(N_ref))= 0;
    if is_alb_opt
        im_target = im;
        im_target(isnan(N_ref))= 0;
    end
    subplot(1,3,1);imshow(im_target)
    title('Target Image')
    subplot(1,3,2);imshow(render_model_noGL(n_ref,l_est_nonamb_lin,alb_ref,0))
    title('Reference')
    subplot(1,3,3);imshow(render_model_noGL(n_new,l_est_nonamb_lin,alb_render,0))
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
    im_diff = alb./N_ref_new.*(l_est_nonamb_lin(2)*p_new+l_est_nonamb_lin(3)*q_new-l_est_nonamb_lin(4))-im*255;
    im_diff(isnan(im_diff))=0;
    im_diff(eye_mask==0) = 0;
    subplot(1,2,2);
    imagesc(abs(im_diff)*100/255);
    title('error in rendered image (Scale of 0-100)')
    
    drawnow
    % depth = estimate_depth(N_ref,alb_ref,im,dmap_ref,l,50);
    
    if mod(count,3)==0 && i<n
        f=figure;
        count = 0;h
        
    end
    count = count+1;
end