function [ z_o,alb_o,L_sh_o ] = intrinsic_decomposition( im,z_ref,labels,mask,rad )
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
[D,S,chrom_p] = seperate_specular(im,mask,labels);
n=2;
nans{1} = isnan(z_ref);
z_ref_iall{1} = z_ref;
sz{1} = size(im);
im(isnan(repmat(z_ref,1,1,3))) = nan;
im_i_all{1} = im;
mask_i_all{1} = mask;
for i=2:n
    z_ref_iall{i}=scale_down_nan(z_ref_iall{i-1})/2;
    nans{i} = isnan(z_ref_iall{i});
    im_i_all{i} = scale_down_nan(im_i_all{i-1});
    sz{i} = size(im_i_all{i});
    mask_i_all{i} = impyramid(mask_i_all{i-1},'reduce');
    
end
for i=n:-1:1
    sc = 2^(-(i-1));
    sz_i = sz{i}(1:2);
    im_i = im_i_all{i};
    if exist('labels_i')
        z_ref_i = scale_up_nan(z_o,sz_i(1:2),nans{i})*2;
        z_ref_i = z_ref_iall{i};
        labels_i = imresize(labels,size(z_ref_i),'nearest');
        mask_i = imresize(mask,size(z_ref_i),'nearest');
        D_i = imresize(D,size(z_ref_i));
        S_i = imresize(S,size(z_ref_i));
    else
        z_ref_i = z_ref_iall{i};
        labels_i = imresize(labels,sc,'nearest');
%         mask_i = imresize(mask,sc,'nearest');
        mask_i = mask_i_all{i};
        D_i = imresize(D,sc);
        S_i = imresize(S,sc);
    end
    tic
    is_face = ~isnan(z_ref_i);
    is_face = remove_bad_boundary(is_face);

    rad_i = rad/sc;
    
    
    
    is_face3 = repmat(is_face,1,1,3);
    
    
    labels3 = repmat(labels_i,1,1,3);
    w_o = [0 0 -1]';
    
    % env map, illumination
    env_path = '../datasets/sIBL jpegs';
    [L_sh_mean, L_sh_var] = compute_sibl_gauss(env_path,8,1);
    L_sh_mean = L_sh_mean/norm(L_sh_mean);
    % alb
    merl_path = '../datasets/merl eth/';
    [alb_mean_cell, alb_var_cell] = load_merl_alb(merl_path);
    [n_mean_cell, n_var_cell,rho_mean_cell,rho_var_cell] = load_merl_shading(merl_path);
    
    
    n_mean = cell2im(n_mean_cell,labels_i,is_face);
    n_var= cell2im(n_var_cell,labels_i,is_face);
    rho_mean = cell2im(rho_mean_cell,labels_i,is_face);
    rho_var= cell2im(rho_var_cell,labels_i,is_face);
    alb_mean = cell2im(alb_mean_cell,labels3,is_face3);
    alb_var = alb_var_cell;
    
    
    
    rho_init = rho_mean;
    delta_a_init = [0.07; 0.23; 0.48];
    delta_s1_init = [1.4; 1.2; 1.3];
    P_init = 1 ;
    
    L_sh_init = L_sh_mean*2;
    alb_init = alb_mean;
    n_init =  n_mean;
    if i<n
        L_sh_init = L_sh_o;
        alb_t = imresize(alb_o,[size(is_face,1),size(is_face,2)]);
        alb_init = alb_t(is_face3);
    end
    [ z_o,alb_o,n_o,rho_o,L_sh_o,chrom_o,delta_a_o,delta_s1_o,P_o ] = ...
        skin_model_optimization( is_face,L_sh_mean,L_sh_var,z_ref_i,chrom_p,...
        alb_init,n_init,rho_init,L_sh_init,D_i,S_i,w_o,labels_i,im_i,alb_mean_cell,...
        alb_var,n_mean,n_var,rho_mean,rho_var,delta_a_init,delta_s1_init,P_init,rad_i,mask_i,1);
    figure;imshow(alb_o.^(1/2.2)*2);
    alb_init = alb_o(is_face3);
    n_init = n_o(is_face);
    rho_init = rho_o(is_face);
    
%     [ z_o,alb_o,n_o,rho_o,L_sh_o,chrom_o,delta_a_o,delta_s1_o,P_o ] = ...
%         skin_model_optimization( is_face,L_sh_mean,L_sh_var,z_o,chrom_o,...
%         alb_init,n_init,rho_init,L_sh_o,D_i,S_i,w_o,labels_i,im_i,alb_mean_cell,alb_var,...
%         n_mean,n_var,rho_mean,rho_var,delta_a_o,delta_s1_o,P_o,rad_i,mask_i,1);
%     figure;imshow(alb_o.^(1/2.2)*2);
    figure;plot(L_sh_o);
    toc
end

%
% subplot(2,1,2);tic
% for i=1:100
% costfn_opti(params);
% end
% toc




end

