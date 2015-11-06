function [ cost,J ] = costfn_skin_model(params, z_ref,chrom_p,D,S,w_o,labels3,wr,wc,alb_mean,alb_var,L_sh_mean,L_sh_var,n_mean,n_var,rho_mean,rho_var,iz_diff,sz,is_face,is_face3,r )
%COSTFN_SKIN_MODEL Summary of this function goes here
%   Detailed explanation goes here
try
z_v = params(1:sz(1));
z2 = nan(size(is_face));
z2(is_face) = z_v;
alb_v = params((sz(1)+1):4*sz(1));
alb = zeros(size(labels3));
alb(is_face3) = alb_v(:);
theta = params(4*sz(1)+1:6*sz(1));
P = params(end);
L_sh = params(6*sz(1)+1:6*sz(1)+sz(2));
L_sh = L_sh*P;
chrom  = params(6*sz(1)+sz(2)+1:6*sz(1)+sz(2)+3);
delta_a = params(end-6:end-4);
delta_s1 = params(end-3:end-1);
n = theta(1:end/2);
rho_s = theta(end/2+1:end);
N = normal_from_depth(z2); 
N_v = reshape(N(is_face3),[],3);
[~,S_r,D_r]=render_DS(N_v',w_o,n,rho_s,alb,r,L_sh,is_face,delta_a,delta_s1,chrom);
%% Data cost
S_r = S_r(is_face);
D_r = D_r(is_face3);
cost_data_spec= S_r(:)-S(:);
cost_data_spec = cost_data_spec*0;
cost_data_diffuse = D_r(:)-D(:);
% cost_data = sum(cost_data_diffuse) + sum(cost_data_spec.^2);
%% Albedo cost
% prior

cost_alb_p = [];
for reg=1:10
    alb_r = reshape(alb(labels3==reg),[],3);
    mean_center = alb_r'-repmat(alb_mean{reg},1,size(alb_r,1));
    cost_alb_p = [cost_alb_p dot(mean_center,(alb_var{reg}\mean_center))];
end
% because later on it will be squared and summed by the optimizer itself.
cost_alb_p = (cost_alb_p').^0.5;

% smoothness
% right neighbour
alb1 = reshape(alb(:,1:end-1,:)-alb(:,2:end,:),[],3);
neg_cost1 = sum( max(-reshape(alb(:,1:end-1,:),[],3),0),2) + sum( max(-reshape(alb(:,2:end,:),[],3),0),2);
alb1 = sum(alb1'.^2).*wc{1}.*wr{1}' + neg_cost1' *10000;
alb1 = [zeros(size(alb,1),1) reshape(alb1,size(alb,1),size(alb,2)-1)];
alb1 = alb1(is_face);
% left neighbour
% alb2 = reshape(alb(:,2:end,:)-alb(:,1:end-1,:),[],3);
% alb2 = sum(alb2'.^2).*wc{1}.*wr{1}';
% alb2 = [reshape(alb2,size(alb,1),size(alb,2)-1) zeros(size(alb,1),1) ];
% alb2 = alb2(is_face);
alb2 = [];
% bot neighbour
alb3 = reshape(alb(1:end-1,:,:)-alb(2:end,:,:),[],3);
neg_cost3 = sum( max(-reshape(alb(1:end-1,:,:),[],3),0),2) + sum( max(-reshape(alb(2:end,:,:),[],3),0),2);
alb3 = sum(alb3'.^2).*wc{2}.*wr{2}' + neg_cost3' *10000;
alb3 = [reshape(alb3,size(alb,1)-1,size(alb,2)); zeros(1,size(alb,2))];
alb3 = alb3(is_face);

% top neighbour
% alb4 = reshape(alb(2:end,:,:)-alb(1:end-1,:,:),[],3);
% alb4 = sum(alb4'.^2).*wc{2}.*wr{2}' ;
% alb4 = [zeros(1,size(alb,2)); reshape(alb4,size(alb,1)-1,size(alb,2))];
% alb4 = alb4(is_face);
alb4 = [];
cost_alb_s = [alb1; alb2; alb3; alb4].^0.5;
% alb_v_r = alb_v(1:end/3);
% alb_v_g = alb_v(end/3+1:2*end/3);
% alb_v_b = alb_v(2*end/3+1:end);
% cost_alb_s = sum(sum((repmat(alb_v(iz_diff(10,:))',4,1)-alb_v(iz_diff([6 9 11 14],:))).^2));


%% Blinn Phong smoothness costs skipped for now
% smoothness
cost_rho_s = repmat(rho_s(iz_diff(10,:))',4,1)-rho_s(iz_diff([6 9 11 14],:));
cost_rho_s = cost_rho_s(:);
cost_n_s = repmat(n(iz_diff(10,:))',2,1)-n(iz_diff([11 14],:));
cost_n_s =cost_n_s(:) ;
% prior
cost_n_p = ((n-n_mean).^2./n_var).^0.5;
%% Geometry cost
% prior
cost_geo_p = z_v-z_ref;
% smoothness
cost_geo_s = repmat(z_v(iz_diff(10,:))',2,1)-z_v(iz_diff([11 14],:));
cost_geo_s =cost_geo_s(:);
%% BSSRDF cost
% absorption parameters
sc = -[0.0305 0.1236 0.1640]';
ic = [0.1 0.331 0.605]';
cost_c = sqrt(sum((delta_a-sc.*delta_s1-ic).^2));
%% Illumination cost
cost_illum_p = chrom-chrom_p;
L_c = L_sh/P-L_sh_mean;
cost_illum_sh = (L_c'/L_sh_var*L_c).^0.5;
cost_illum_n = (norm(chrom)-1) + (norm(L_sh/P)-1) + sum(sum(max(-sh2sph(L_sh/P,[50 50]),0)))/5*0 ;

lamb_fp = sqrt(0);
lamb_hs = sqrt(0);
lamb_zs = sqrt(0);
lamb_ln = sqrt(0);
lamb_fs = sqrt(0);
lamb_lp = sqrt(0);
lamb_lsh = sqrt(0);
lamb_hp = sqrt(0);
lamb_hc = sqrt(0);
lamb_zp = sqrt(0);

% lamb_fp = sqrt(10);
% lamb_fp = sqrt(0.0001);
% lamb_hs = sqrt(1);
% lamb_zs = sqrt(0.01);
% lamb_ln = sqrt(10);
% lamb_fs = sqrt(10);
% lamb_lp = sqrt(5);
% lamb_lsh = sqrt(0.00001);
% lamb_hp = sqrt(0.5);
% lamb_hc = sqrt(10);
% lamb_zp = sqrt(0.0001);
lamb_fp = sqrt(0.001);
lamb_hs = sqrt(1000);
lamb_zs = sqrt(0);
lamb_ln = sqrt(10);
lamb_fs = sqrt(1000);
lamb_lp = sqrt(0.01);
lamb_lsh = sqrt(1);
lamb_hp = sqrt(0.5);
lamb_hc = sqrt(0.1);
lamb_zp = sqrt(0.01);
cost = [cost_data_diffuse*20; cost_data_spec ;cost_alb_p*lamb_fp;...
    cost_alb_s*lamb_fs; cost_n_p*lamb_hp; cost_n_s*lamb_hs;...
    cost_c*lamb_hc; cost_geo_p*lamb_zp; cost_geo_s*lamb_zs; ...
    cost_illum_p*lamb_lp; cost_illum_sh*lamb_lsh*0; cost_illum_n*lamb_ln];

if sum(~isreal(cost))>0
    disp(1)
end
catch err
    disp(1)
end
% cost = cost*0;
% cost(27495)=1-params(13748);
end

