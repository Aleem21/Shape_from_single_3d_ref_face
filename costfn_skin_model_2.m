function [ cost ] = costfn_skin_model_2(params,D,S,w_o,labels3,sz,is_face,is_face3,r)
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

cost = [cost_data_diffuse*20; cost_data_spec];
if sum(~isreal(cost))>0
    disp(1)
end
catch err
    error('Error here')
end

% cost = cost*0;
% cost(27495)=1-params(13748);
end

