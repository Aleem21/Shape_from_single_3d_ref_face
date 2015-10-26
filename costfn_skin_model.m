function [ cost ] = costfn_skin_model(params, z_ref,chrom_p,D,S,w_o,labels3,wr,wc,alb_mean,alb_var,L_sh_mean,L_sh_var,n_mean,n_var,rho_mean,rho_var,iz_diff,sz,is_face,is_face3,r )
%COSTFN_SKIN_MODEL Summary of this function goes here
%   Detailed explanation goes here
z_v = params(1:sz(1));
z2 = nan(size(is_face));
z2(is_face) = z_v;
alb_v = params(sz(1)+1:4*sz(1));
alb = zeros(size(labels3));
alb(is_face3) = alb_v(:);
theta = params(4*sz(1)+1:6*sz(1));
L_sh = params(6*sz(1)+1:6*sz(1)+sz(2));
L_sh = L_sh./norm(L_sh);
chrom  = params(6*sz(1)+sz(2)+1:6*sz(1)+sz(2)+3);
delta_a = params(end-5:end-3);
delta_s1 = params(end-2:end);
n = theta(1:end/2);
rho_s = theta(end/2+1:end);
N = normal_from_depth(z2); 
N_v = reshape(N(is_face3),[],3);
[~,S_r,D_r]=render_DS(N_v',w_o,n,rho_s,alb,r,L_sh,is_face,delta_a,delta_s1);
%% Data cost
S_r = S_r(is_face);
D_r = D_r(is_face3);
cost_data_diffuse = S_r(:)-S(:);
cost_data_spec = D_r(:)-D(:);
cost_data = sum(cost_data_diffuse) + sum(cost_data_spec.^2);
%% Albedo cost
% prior
cost_alb_p = 0;
for reg=1:10
    alb_r = reshape(alb(labels3==reg),[],3);
    mean_center = alb_r'-repmat(alb_mean{reg},1,size(alb_r,1));
    cost_alb_p = cost_alb_p + sum(dot(mean_center,(alb_var{reg}\mean_center)));
end

% smoothness
% right neighbour
alb1 = reshape(alb(:,1:end-1,:)-alb(:,2:end,:),[],3);
alb1 = sum(alb1'.^2).*wc{1}.*wr{1}';
% left neighbour
alb2 = reshape(alb(:,2:end,:)-alb(:,1:end-1,:),[],3);
alb2 = sum(alb2'.^2).*wc{1}.*wr{1}';
% bot neighbour
alb3 = reshape(alb(1:end-1,:,:)-alb(2:end,:,:),[],3);
alb3 = sum(alb3'.^2).*wc{2}.*wr{2}';
% top neighbour
alb4 = reshape(alb(2:end,:,:)-alb(1:end-1,:,:),[],3);
alb4 = sum(alb4'.^2).*wc{2}.*wr{2}';

cost_alb_s = sum(alb1)+sum(alb2)+sum(alb3)+sum(alb4);
% alb_v_r = alb_v(1:end/3);
% alb_v_g = alb_v(end/3+1:2*end/3);
% alb_v_b = alb_v(2*end/3+1:end);
% cost_alb_s = sum(sum((repmat(alb_v(iz_diff(10,:))',4,1)-alb_v(iz_diff([6 9 11 14],:))).^2));


%% Blinn Phong smoothness costs skipped for now
% smoothness
cost_rho_s = sum(sum((repmat(rho_s(iz_diff(10,:))',4,1)-rho_s(iz_diff([6 9 11 14],:))).^2));
cost_n_s = sum(sum((repmat(n(iz_diff(10,:))',4,1)-n(iz_diff([6 9 11 14],:))).^2));

% prior
cost_n_p = (n-n_mean).^2./n_var;
%% Geometry cost
% prior
cost_geo_p = sum((z_v-z_ref).^2);
% smoothness
cost_geo_s = sum(sum((repmat(z_v(iz_diff(10,:))',4,1)-z_v(iz_diff([6 9 11 14],:))).^2));

%% BSSRDF cost
n_mat = zeros(size(z2));
n_mat(is_face) = n;
rho_mat = zeros(size(z2));
rho_mat(is_face) = rho_s;
% prior
cost_n_p = (n-n_mean).^2./n_var;

% absorption parameters
sc = rand(3,1);
ic = rand(3,1);
cost_c = sum((delta_a-sc.*delta_s1-ic).^2);
%% Illumination cost
cost_illum_p = sum((chrom-chrom_p).^2);
cost_illum_sh = (L_sh-L_sh_mean)'*(L_sh_var\(L_sh-L_sh_mean));
cost_illum_n = (norm(chrom)-1)^2;


lamb_fp = 10;
lamb_hs = 10;
lamb_zs = 10;
lamb_ln = 10;
lamb_fs = 1;
lamb_lp = 1;
lamb_lsh = 1;
lamb_hp = 0.5;
lamb_hc = 0.1;
lamb_zp = 0.01;

cost = cost_data + cost_alb_p * lamb_fp+ cost_alb_s * lamb_fs +...
    cost_geo_p * lamb_zp + cost_geo_s * lamb_zs + cost_illum_p * lamb_lp +...
    cost_illum_sh * lamb_lsh + cost_illum_n * lamb_ln;
cost = [cost_data_diffuse; cost_data_spec;alb1'; alb2'; alb3'; alb4'; cost_alb_p; ];
if nargout >1
    r = 1:numel(cost_data_diffuse);
%     c = 
end
end

