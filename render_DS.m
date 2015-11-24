function [ I,S,I_D,alpha,Pdr_filt,Yi,Ft_o,E_lm_out,d_S_n,l,Y_N ] = render_DS( N_v,w_o,n,rho_s,alb,r,Light_sh,is_face,delta_a,delta_s1,chrom )
%RENDER_DS Summary of this function goes here
%   Detailed explanation goes here

%% Definitions
% w_o                       - viewing direction
% rho_s                     - specular scale
% n                         - specular exponent
% surface normals

%% SH approach
% N_wi_theta = acos(N_dot_wi);
N_wo_theta = acos(dot(N_v,repmat(w_o,1,size(N_v,2)) ));

aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

[rec_pts_az,rec_pts_el] = cart2sph(N_v(1,:),N_v(2,:),N_v(3,:));
rec_pts = aziElev2aziIncl([rec_pts_az' rec_pts_el']);

R_v = 2*repmat(dot(repmat(w_o,1,size(N_v,2)),N_v),3,1).*N_v-repmat(w_o,1,size(N_v,2));
[rot_pts_az,rot_pts_el] = cart2sph(R_v(1,:),R_v(2,:),R_v(3,:));
rot_pts = aziElev2aziIncl([rot_pts_az' rot_pts_el']);

Y_N = getSH(8, [rec_pts; rot_pts], 'real');

n_air = 1;
n_skin = 1.38;
% Diffuse (BSSRDF)
Nord = 8;
[~, tdesign_grid] = getTdesign(2*Nord);
tdesign_grid = aziElev2aziIncl(tdesign_grid); % convert to azi-incl
% Ft_sh = max(cos(tdesign_grid(:,2)),0).*fresnel_trans(tdesign_grid(:,2),n_air,n_skin,1);
% Ft_sh = cos(tdesign_grid(:,2))>0;
Ft_sh = (cos(tdesign_grid(:,2))>0).*fresnel_trans(tdesign_grid(:,2),n_air,n_skin,1);

Ft_sh_coeff = directSHT(Nord, Ft_sh, tdesign_grid, 'real', []);
l = zeros((Nord+1)^2,1);
for i=0:Nord
    l(i^2+1:(i+1)^2) = i;
end
const = sqrt(4*pi./(2*l+1));

E_lm = const.*Ft_sh_coeff;
for i=0:Nord
    E_lm(i^2+1:(i+1)^2) = E_lm(median([i^2+1 (i+1)^2]));
end

if numel(Light_sh)<numel(E_lm)
    Light_sh(numel(E_lm)) = 0;
end
in = inverseSHT(E_lm.*Light_sh, rec_pts, 'real');

Ft_o_v = fresnel_trans(N_wo_theta,n_skin,n_air,2);
in_integ = zeros(size(is_face));
in_integ(is_face) = in;
in_integ(isnan(in_integ))=0;
Ft_o_v(isnan(Ft_o_v))=0;
Ft_o = zeros(size(is_face));
Ft_o(is_face) = Ft_o_v;
for clr=1:3
    [~,Pdr_filt{clr}]=bssdf([],[],r,1,max(delta_a,0),max(delta_s1,0));
    Pdr_filt{clr} = Pdr_filt{clr}(2:end-1,2:end-1);
%     D(:,:,clr) = 4*conv2(Ft_o.*in_integ.*alb(:,:,clr),Pdr_filt{clr}*r^2,'same')*chrom(clr,:);
    D(:,:,clr) = 4*conv2(in_integ,Pdr_filt{clr}*r^2,'same')*chrom(clr,:).*alb(:,:,clr).*Ft_o;
end
if nargout > 3
    alpha = in_integ.*Ft_o.*r^2;
    alpha = alpha(is_face);
    Yi = Y_N(1:size(rec_pts,1),:);
    E_lm_out = E_lm;
end
D(D<0) = 0;

%% SH specular
Nord = 8;
l = zeros((Nord+1)^2,1);
for i=0:Nord
    l(i^2+1:(i+1)^2) = i;
end
l = repmat(l,1,numel(n));
n2 = repmat(n(:)',size(l,1),1);
A_lm = exp(-l.^2/2./n2);

E_lm = repmat(Light_sh,1,size(A_lm,2)).*A_lm;
% S_v = rho_s'.*inverseSHT_my(E_lm, rot_pts, 'real',Y_N(size(rec_pts,1)+1:end,:));
Y_N = getSH(Nord,rot_pts,'real');
% S_v = rho_s'.*inverseSHT_my(E_lm, rot_pts, 'real');

S_v = rho_s'.*inverseSHT_my(E_lm, [], 'real',Y_N);
S = zeros(size(is_face));
S(is_face) = S_v;
S(S<0) = 0;
if nargout >3
    E_lm_d_S_n = E_lm.*l.^2./(2*n2.^2);
    d_S_n = rho_s'.*inverseSHT_my(E_lm_d_S_n, rot_pts, 'real');
end
I = repmat(S,1,1,3) + D;

I_D = D;
end

