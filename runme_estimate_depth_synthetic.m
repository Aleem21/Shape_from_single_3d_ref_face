clear variables
% close all
import plot3D_helper.label_axis

talk = 0;
%% set initial variables
% folder_path = '.\data\l_test\';

%% generate ref depth map
center = 1;
if center
    s = 150;
    cRes = 200 ;rRes = 200;
    Rpose = makehgtform('scale',[rRes/s cRes/s  1/s]);
    Scale = [rRes cRes];
else
    s = 150;
    Rpose = makehgtform('scale',1/s);
    cRes = 200 ;rRes = 200;
    Scale = [rRes cRes]/2;
end

ply_path = '.\data\sphere.ply';
% Rpose = makehgtform('scale',[1 1  1/150]);
% % Rpose = makehgtform('scale',[1/150]);
% cRes = 200 ;rRes = 200;
% Scale = [150 150];
% % Scale = [150 150];
% center = 1;
% center = 0;
Rpose = eye(4);
ply_path = '.\data\ref_model.ply';
center = 0;
% [dmap_ref, n_ref, N_ref] = generate_ref_depthmap(ply_path,Scale, talk, cRes, rRes, Rpose,[],center);
[dmap_ref, n_ref, N_ref] = generate_ref_depthmap_synth(100,1);
figure; surf(dmap_ref,'edgealpha',0.5);title('Original Model')


%% generate ref albedo
alb_ref = ~isnan(N_ref);

%% render image
l_ren = [0;-0.6;0.5;0.8]/2;
% im = render_model(ply_path,l_ren,talk,Rpose, 150,150);
im = render_model_noGL( n_ref, l_ren, alb_ref,1);
title('image')
%% estimate lighting
talk = 0;
l_est = estimate_lighting(n_ref, alb_ref, im,4);
% c4 = render_model('./data/sphere.ply',l_est,talk,Rpose, 150,150);
c4 = render_model_noGL(n_ref,l_est,alb_ref,talk);
noise = randn(size(dmap_ref))/8;
depth = estimate_depth(N_ref,alb_ref,im,dmap_ref+noise,l_est,0.01);

figure; surf(depth,'edgealpha',0.5);title('reconstructed')
figure; surf(depth-dmap_ref,'edgealpha',0.5);title('error')
error = (depth-dmap_ref).^2;
error(isnan(error)) = 0;
signal_energy = dmap_ref.^2;
signal_energy(isnan(signal_energy)) = 0;

fprintf('Percentage squared error per pixel = %d %%\n',sum(error(:))/sum(signal_energy(:))*100);
fprintf('Singal energy = %d %%\n',sum(signal_energy(:)));
fprintf('Noise energy = %d %%\n',sum(sum(noise.^2)));