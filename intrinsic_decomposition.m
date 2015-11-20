function [ im_o ] = intrinsic_decomposition( im,z_ref,labels,mask,is_face,rad )
%INTRINSIC_DECOMPOSITION Summary of this function goes here
%   Detailed explanation goes here
% complex model
rho_s = zeros(size(im,1),size(im,2))+1;
expo = zeros(size(im,1),size(im,2))+300;
% 
i_path = '../datasets/sIBL/Sierra_Madre_B/Sierra_Madre_B_HiRes_TMap.jpg';

% i_env = imread(i_path);
% i_env = rgb2gray(i_env);

% sh=env2sh(i_env,10,0);
% sh = sh/norm(sh);
% [I,S,D]=render_DS(normal_from_depth(z_ref),[0;0;-1],expo,rho_s,im*0+1,rad,sh);
is_face3 = repmat(is_face,1,1,3);

[D,S,chrom_p] = seperate_specular(im,mask,labels);

labels3 = repmat(labels,1,1,3);
w_o = [0 0 -1]';

% env map, illumination
env_path = '../datasets/sIBL jpegs';
[L_sh_mean, L_sh_var] = compute_sibl_gauss(env_path,8,1);
L_sh_mean = L_sh_mean/norm(L_sh_mean);
% alb
merl_path = '../datasets/merl eth/';
[alb_mean_cell, alb_var_cell] = load_merl_alb(merl_path);
[n_mean_cell, n_var_cell,rho_mean_cell,rho_var_cell] = load_merl_shading(merl_path);


n_mean = cell2im(n_mean_cell,labels,is_face);
n_var= cell2im(n_var_cell,labels,is_face);
rho_mean = cell2im(rho_mean_cell,labels,is_face);
rho_var= cell2im(rho_var_cell,labels,is_face);
alb_mean = cell2im(alb_mean_cell,labels3,is_face3);
alb_var = alb_var_cell;



alb_init = alb_mean;
n_init = n_mean;
rho_init = rho_mean;
delta_a_init = [0.07; 0.23; 0.48];
delta_s1_init = [1.4; 1.2; 1.3];
P_init = 1 ;


L_sh_init = L_sh_mean;
[ z_o,alb_o,n_o,rho_o,L_sh_o,chrom_o,delta_a_o,delta_s1_o,P_o ] = ...
    skin_model_optimization( is_face,L_sh_mean,L_sh_var,z_ref,chrom_p,...
    alb_init,n_init,rho_init,L_sh_init,D,S,w_o,labels,im,alb_mean_cell,...
    alb_var,n_mean,n_var,rho_mean,rho_var,delta_a_init,delta_s1_init,P_init,rad,0);
figure;imshow(alb_o*P_o);
alb_init = alb_o(is_face3);
n_init = n_o(is_face);
rho_init = rho_o(is_face);

[ z_o,alb_o,n_o,rho_o,L_sh_o,chrom_o,delta_a_o,delta_s1_o,P_o ] = ...
    skin_model_optimization( is_face,L_sh_mean,L_sh_var,z_o,chrom_o,...
    alb_init,n_init,rho_init,L_sh_o,D,S,w_o,labels,im,alb_mean_cell,alb_var,...
    n_mean,n_var,rho_mean,rho_var,delta_a_o,delta_s1_o,P_o,rad,1);
figure;imshow(alb_o*50);


subplot(2,1,2);tic
for i=1:100
costfn_opti(params);
end
toc




end

