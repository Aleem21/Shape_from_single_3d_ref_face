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
delta_a = [0.07; 0.23; 0.48];
delta_s1 = [1.4; 1.2; 1.3];
L_sh_init = L_sh_mean;
P = 32 ;
sz = [sum(is_face(:)) numel(L_sh_init)];
params = [z_init(:); alb_init; theta_init(:); L_sh_init(:); chrom_init(:); delta_a; delta_s1; P];
labels3 = labels3.*is_face3;


% Jacobian pattern
sz_fc = sum(is_face(:));
[r_face,c_face] = find(is_face);
[bboxc, bboxr] = meshgrid(-5:5,-5:5);
[boxc, boxr] = meshgrid(-2:1,-2:1);

[boxc3, boxr3] = meshgrid(-1:1,-1:1);
sub2ind_face = zeros(size(is_face));
sub2ind_face(is_face) = 1:sum(is_face(:));
iz_diff = zeros(15,sum(is_face(:)));
iz_spec = zeros(3,sum(is_face(:)));
iz_alb = zeros(9,sum(is_face(:)));
for i=1:sum(is_face(:))
%     belems4x4 = sub2ind_face(sub2ind(size(is_face),min(max(bboxr(2:end)+r_face(i),1),size(is_face,1)),min(max(bboxc(2:end)+c_face(i),1),size(is_face,2))));
    elems4x4 = sub2ind_face(sub2ind(size(is_face),min(max(boxr(2:end)+r_face(i),1),size(is_face,1)),min(max(boxc(2:end)+c_face(i),1),size(is_face,2))));

    elems3x3 = sub2ind_face(sub2ind(size(is_face),boxr3(:)+r_face(i),boxc3(:)+c_face(i)));
    iz_spec(:,i) = elems4x4([6 9 10]) + elems4x4([14 11 10]).*(elems4x4([6 9 10])==0);
    elems4x4(elems4x4==0) = repmat(elems4x4(10),sum(elems4x4==0),1);
   
    iz_diff(:,i) = elems4x4;
    elems3x3(elems3x3==0) = repmat(elems3x3(5),sum(elems3x3==0),1);
    iz_alb(:,i) = elems3x3;
end

% Diffuse data term
P_offset = numel(params);
chrom_offset = 6*sz_fc  + numel(L_sh_mean);
sh_offset = 6*sz_fc;

alb_offset = sz_fc;
r_1_1 = repmat(1:size(iz_diff,2),15+9+1+3+numel(L_sh_init),1);
c_1_1 = [iz_diff; alb_offset + iz_alb ; ones(1,size(iz_diff,2))*P_offset;repmat((1:3)',1,size(iz_diff,2))+chrom_offset; repmat((1:numel(L_sh_init))',1,size(iz_diff,2))+sh_offset; ];
r_1_2 = r_1_1(end) + repmat(1:size(iz_diff,2),15+9+1+3+numel(L_sh_init),1);
c_1_2 = [iz_diff; alb_offset + iz_alb + sz_fc ; ones(1,size(iz_diff,2))*P_offset;repmat((1:3)',1,size(iz_diff,2))+chrom_offset; repmat((1:numel(L_sh_init))',1,size(iz_diff,2))+sh_offset;];
r_1_3 = r_1_2(end) + repmat(1:size(iz_diff,2),15+9+1+3+numel(L_sh_init),1);
c_1_3 = [iz_diff; alb_offset + iz_alb + sz_fc*2 ; ones(1,size(iz_diff,2))*P_offset;repmat((1:3)',1,size(iz_diff,2))+chrom_offset; repmat((1:numel(L_sh_init))',1,size(iz_diff,2))+sh_offset;];

% specular data term
r_2 = repmat(1:size(iz_spec,2),3+1+1+numel(L_sh_mean),1) + r_1_3(end);
% % r_2 = r_2*0 + r_2(end);
c_2 = [iz_spec(:); iz_spec(3,:)'+sz_fc*4; iz_spec(3,:)'+sz_fc*5; repmat((1:numel(L_sh_mean))'+sz_fc*6,sz_fc,1)];
% % c_r = c_2*0 + c_2(end);
% r_2=[;
% albedo prior
r_3 = r_2(end) + repmat(1:size(iz_spec,2),3,1);
% r_3 = r_1_3(end) + repmat(1:size(iz_spec,2),3,1);
c_3 = alb_offset + [iz_diff(10,:); iz_diff(10,:)+sz_fc; iz_diff(10,:)+2*sz_fc];


% Albedo smoothing term
%right
r_4 = r_3(end)+ repmat(1:sz_fc,2,1);
c_4 = alb_offset + iz_diff([10 14],:);
%left
% r_5 = r_4(end)+ repmat(1:sz_fc,2,1);
% c_5 = alb_offset + iz_diff([10 6],:);
r_5 = []; c_5 = [];
%bot
r_6 = r_4(end)+ repmat(1:sz_fc,2,1);
c_6 = alb_offset + iz_diff([11 6],:);
%top
% r_7 = r_6(end)+ repmat(1:sz_fc,2,1);
% c_7 = alb_offset + iz_diff([9 6],:);
r_7  =[]; c_7 = [];
% n prior term
n_offset = 4*sz_fc;
r_8 = r_6(end)+ (1:sz_fc);
c_8 = n_offset + (1:sz_fc);

% n smooting terms
r_9 = r_8(end) + repmat(1:2*sz_fc,2,1);
c_9 = n_offset + iz_diff([10 11 10 14],:);
r_9 = r_9(:);
c_9 = c_9(:);
% absorption term
abs_offset = 6*sz_fc + numel(L_sh_mean) + 3;
r_10 = r_9(end) + ones(1,6);
c_10 = abs_offset + (1:6);

% geometry prior
geo_offset = 0;
r_11 = r_10(end)+ (1:sz_fc);
c_11 = geo_offset + (1:sz_fc);

% geometry smoothing
r_12 = r_11(end) + repmat(1:2*sz_fc,2,1);
c_12 = geo_offset + iz_diff([10 11 10 14],:);

% chrom prior
r_13 = r_12(end) + (1:3);
c_13 = chrom_offset + (1:3);

% SH prior
r_14 = r_13(end) + ones(1,numel(L_sh_mean));
c_14 = sh_offset + (1:numel(L_sh_mean));

% chrom norm
r_15 = r_14(end) + ones(1,3+numel(L_sh_init));
c_15 = [chrom_offset + (1:3) sh_offset + (1:numel(L_sh_init))];
n = r_15(end);
m = numel(params);
J = sparse([r_1_1(:);r_1_2(:);r_1_3(:);r_2(:);r_3(:);r_4(:);r_5(:);r_6(:)...
    ;r_7(:);r_8(:);r_9(:);r_10(:);r_11(:);r_12(:);r_13(:);r_14(:);r_15(:)],...
    [c_1_1(:);c_1_2(:);c_1_3(:);c_2(:);c_3(:);c_4(:);c_5(:);c_6(:)...
    ;c_7(:);c_8(:);c_9(:);c_10(:);c_11(:);c_12(:);c_13(:);c_14(:);c_15(:)],...
    1,n,m);
J(:,1:sz_fc) = 0;
J(J>0) = 1;
% J(:,1:sz_fc)=0;
costfn_opti = @(params)costfn_skin_model(params, z_ref(is_face),chrom_p,D(is_face3),S(is_face),w_o,...
    labels3,wr,wc,alb_mean,alb_var,L_sh_mean,L_sh_var,n_mean,n_var,rho_var,rho_mean,iz_diff,sz,is_face,is_face3,rad );
% J = sparse(27495,13747,-1,numel(cost),numel(params));
algo = 'levenberg-marquardt';
% algo =  {'levenberg-marquardt',.005};
algo = 'trust-region-reflective';
options = optimset('MaxIter',50,'JacobPattern',J,'Display','iter-detailed','Algorithm',algo);
tic
optimum = lsqnonlin(costfn_opti,params,[],[],options);
toc
figure;plot(params-optimum);
disp(1);
z_o = nan(size(z_ref));
z_o(is_face) = optimum(1:sz_fc);
alb_o = zeros(size(im));
alb_o(is_face3) = optimum(sz_fc+1:4*sz_fc);
n_o = nan(size(z_ref));
n_o(is_face) = optimum(4*sz_fc+1:5*sz_fc);
rho_o = nan(size(z_ref));
rho_o(is_face) = optimum(5*sz_fc+1:6*sz_fc);
L_sh_o = optimum(6*sz_fc+[1:numel(L_sh_mean)]);
chrom_o = optimum(6*sz_fc+numel(L_sh_mean)+[1:3]);
delta_a_o = optimum(end-6:end-4);
delta_s1_o = optimum(end-3:end-1);
P_o = optimum(end);
% params = [z_init(:); alb_init; theta_init(:); L_sh_init(:); chrom_init(:); delta_a; delta_s1];

% options = optimset('MaxIter',2);
figure;
subplot(2,1,1);
ci = costfn_opti(params).^2;
plot(1:r_1_3(end),ci(1:r_1_3(end)),r_2(1):r_2(end),ci(r_2(1):r_2(end)),...
    r_3(1):r_3(end),ci(r_3(1):r_3(end)),r_4(1):r_8(end),ci(r_4(1):r_8(end)),...
    r_9(1):r_9(end),ci(r_9(1):r_9(end)),r_9(end):r_11(1),ci(r_9(end):r_11(1)),r_11(1):r_15(end),ci(r_11(1):r_15(end)));
subplot(2,1,2);
co = costfn_opti(optimum).^2;
plot(1:r_1_3(end),co(1:r_1_3(end)),r_2(1):r_2(end),co(r_2(1):r_2(end)),...
    r_3(1):r_3(end),co(r_3(1):r_3(end)),r_4(1):r_8(end),co(r_4(1):r_8(end)),...
    r_9(1):r_9(end),co(r_9(1):r_9(end)),r_9(end):r_11(1),co(r_9(end):r_11(1)),r_11(1):r_15(end),co(r_11(1):r_15(end)));
subplot(2,1,2);tic
for i=1:100
costfn_opti(params);
end
toc




end

