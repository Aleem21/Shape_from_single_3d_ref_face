function [ output_args ] = intrinsic_decomposition( im,z_ref,labels,mask,is_face,rad )
%INTRINSIC_DECOMPOSITION Summary of this function goes here
%   Detailed explanation goes here
% complex model
figure;
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
[alb_mean, alb_var] = load_merl_alb(merl_path);
[n_mean1, n_var1,rho_mean1,rho_var1] = load_merl_shading(merl_path);
alb_init = zeros(size(im));
n_mean = zeros(size(im,1),size(im,2));
rho_mean = n_mean;
n_var = n_mean;
rho_var = rho_mean;
for r=1:10
    temp = repmat(alb_mean{r}',sum(labels3(:)==r)/3,1);
    alb_init(labels3==r) = temp(:);
    temp = repmat(n_mean1{r}',sum(labels3(:)==r)/3,1);
    n_mean(labels==r) = temp(:);
    n_var(labels==r) = n_var1{r};
    temp = repmat(rho_mean1{r}',sum(labels3(:)==r)/3,1);
    rho_mean(labels==r) = temp(:);
    rho_var(labels==r) = rho_var1{r};
end
n_var = n_var(is_face);
rho_var = rho_var(is_face);
alb_init = alb_init.*repmat(is_face,1,1,3);
alb_init = alb_init(is_face3);

n_mean = n_mean.*is_face;
n_mean = n_mean(is_face);

rho_mean = rho_mean.*is_face;
rho_mean = rho_mean(is_face);

% wr,wc weights albedo smoothness
wr{1} = (1 + (labels(:,1:end-1)~=labels(:,2:end))*0.5).*is_face(:,1:end-1).*is_face(:,2:end);
wr{1} = wr{1}(:);
wr{2} = (1 + (labels(2:end,:)~=labels(1:end-1,:))*0.5).*is_face(1:end-1,:).*is_face(2:end,:);
wr{2} = wr{2}(:);

chrom = reshape(im,[],3);
chrom = chrom./repmat(sqrt(sum(chrom.^2,2)),1,3);
chrom = reshape(chrom,size(im));
chroml = reshape(chrom(:,1:end-1,:),[],3);
chromr = reshape(chrom(:,2:end,:),[],3);
wc{1} = 1 - (dot(chroml',chromr')<0.9)*0.6;
chromt = reshape(chrom(1:end-1,:,:),[],3);
chromb = reshape(chrom(2:end,:,:),[],3);
wc{2} = 1 - (dot(chromt',chromb')<0.9)*0.6;


theta_init = [n_mean; rho_mean];
chrom_init = chrom_p;
z_init = z_ref(is_face);
% alb_init = 
delta_a = [0.7; 0.23; 0.48];
delta_s1 = [1.4; 1.2; 1.3];
L_sh_init = L_sh_mean;
sz = [sum(is_face(:)) numel(L_sh_init)];
params = [z_init(:); alb_init; theta_init(:); L_sh_init(:); chrom_init(:); delta_a; delta_s1];
labels3 = labels3.*is_face3;

costfn_opti = @(params)costfn_skin_model(params, z_ref(is_face),chrom_p,D(is_face3),S(is_face),w_o,...
    labels3,wr,wc,alb_mean,alb_var,L_sh_mean,L_sh_var,n_mean,n_var,rho_var,rho_mean,sz,is_face,is_face3,rad );
options = optimset('MaxIter',2);
fminunc(costfn_opti,params,options);



% regularization term
% in_inface = (in_face-b_in_full);

% face_inds = find(is_face);
[r_face,c_face] = find(is_face);

[boxc, boxr] = meshgrid(-2:1,-2:1);

sub2ind_face = zeros(size(is_face));
sub2ind_face(is_face) = 1:sum(is_face(:));

iz_reg = zeros(sum(is_face(:)),numel(boxc));
for i=1:sum(is_face(:))
    elems4x4 = sub2ind_face(sub2ind(size(is_face),boxr(:)+r_face(i),boxc(:)+c_face(i)));
    ind_out_of_bounds = find(elems4x4==0);
    iz_reg(i,:) = elems4x4;
end





end

