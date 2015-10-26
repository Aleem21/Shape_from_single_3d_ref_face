function [ I,S,I_D ] = render_DS( N_v,w_o,n,rho_s,alb,r,Light_sh,is_face,delta_a,delta_s1 )
%RENDER_DS Summary of this function goes here
%   Detailed explanation goes here

%% Definitions
% w_o                       - viewing direction
% rho_s                     - specular scale
% n                         - specular exponent
% surface normals
%% Specular part
% % light direction
% L_d = L/norm(L);
% L_mag = norm(L);
% % half angle vector
% H = -(L-w_o)/norm(L_d+w_o);

% N_v(1,:) = reshape(N(:,:,1),1,[]);
% N_v(2,:) = reshape(N(:,:,2),1,[]);
% N_v(3,:) = reshape(N(:,:,3),1,[]);
% H_v = repmat(H,1,size(N_v,2));
% cos_delta = dot(N_v,H_v);
% cos_delta = max(reshape(cos_delta,size(N,1),size(N,2)),0);
% 
% % S = zeros(size(N,1),size(N,2));
% % for i=1:size(S,1)
% %     for j=1:size(S,2)
% %        S(i,j) = rho_s*(n+2)/(2*pi) * 
% %     end
% % end
% N_dot_wi = dot(N_v,repmat(-L_d,1,size(N_v,2)));
% N_dot_wi = max(reshape(N_dot_wi,size(N,1),size(N,2)),0);
% S = rho_s * (n+2)/(2*pi) * cos_delta.^n * L_mag .* N_dot_wi;


%% SH approach
% N_wi_theta = acos(N_dot_wi);
N_wo_theta = acos(dot(N_v,repmat(w_o,1,size(N_v,2)) ));

aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

[rec_pts_az,rec_pts_el] = cart2sph(N_v(1,:),N_v(2,:),N_v(3,:));
rec_pts = aziElev2aziIncl([rec_pts_az' rec_pts_el']);

n_air = 1;
n_skin = 1.38;
% Diffuse (BSSRDF)
Nord = 10;
[~, tdesign_grid] = getTdesign(2*Nord);
tdesign_grid = aziElev2aziIncl(tdesign_grid); % convert to azi-incl
Ft_sh = max(cos(tdesign_grid(:,2)),0).*fresnel_trans(tdesign_grid(:,2),n_air,n_skin,1);
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

% nx = N(:,:,1);
% nx = nx(:)';        %change to row vector
% ny = N(:,:,2);
% ny = ny(:)';
% nz = N(:,:,3);
% nz = nz(:)';
% 
% Y = [ones(size(nx));
%     nx;
%     ny;
%     nz;
%     nx.*ny;
%     nx.*nz;
%     ny.*nz;
%     (nx.^2 - ny.^2);
%     (3*nz.^2-1) ];
if numel(Light_sh)<numel(E_lm)
    Light_sh(numel(E_lm)) = 0;
end
in = inverseSHT(E_lm.*Light_sh, rec_pts, 'real');
% im = (E_lm'.*Light_sh)*Y;
% in = reshape(in,size(N,1),size(N,2));
% figure;
% imshow(in);

% 
Ft_o_v = fresnel_trans(N_wo_theta,n_skin,n_air,2);
in_integ = zeros(size(is_face));
in_integ(is_face) = in;
in_integ(isnan(in_integ))=0;
Ft_o_v(isnan(Ft_o_v))=0;
Ft_o = zeros(size(is_face));
Ft_o(is_face) = Ft_o_v;
for clr=1:3
    [~,Pdr_filt]=bssdf([],[],r,clr,delta_a,delta_s1);
    D(:,:,clr) = 4*conv2(Ft_o.*in_integ,Pdr_filt*r^2,'same');
end

% Specular (Blinn-Phong)
% Nord = 10;
% aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];
% [~, tdesign_grid] = getTdesign(2*Nord);
% tdesign_grid = aziElev2aziIncl(tdesign_grid); % convert to azi-incl
% l = zeros((Nord+1)^2,1);
% for i=0:Nord
%     l(i^2+1:(i+1)^2) = i;
% end
% const = sqrt(4*pi./(2*l+1));
% 
% 
% 
% grid_sph = aziElev2aziIncl(tdesign_grid);
% % [x,y,z] = sph2cart(grid_sph(:,1),grid_sph(:,2),ones(size(grid_sph,1),1));
% % grid_cart = [x';y';z'];
% % H = grid_cart + repmat(w_o,1,size(grid_cart,2));
% % H = H./repmat(sqrt(sum(H.^2)),3,1);
% % generate sh coeff for second term in the integral of blinn-phong. Do so
% % at each point of our grid
% Nord = 10;
% n=1;
% for i=1:size(grid_sph,1)
%     [Nix,Niy,Niz] = sph2cart(grid_sph(i,1),grid_sph(i,2),1);
%     Ni = [Nix;Niy;Niz];
%     sh_c(:,i) = find_blinnp_sh_coeff(Ni,w_o,tdesign_grid,n,Nord);
% end
% % fit spherical harmonics in each sh coefficients, thus we can generate
% % those SH_coefficients from these sets given any normal
% for i=1:size(sh_c,1)
%     sh_c2(:,i) = directSHT(Nord, sh_c(i,:)', tdesign_grid, 'real', []);
% end
% 
% %generate the SH_coefficients of the second term given any normal, from
% %their second level SH's
% sh_c_recov = zeros(size(sh_c2,2),1);
% N_recovery = rand(3,1); N_recovery = N_recovery/norm(N_recovery);
% [rec_az,rec_el] = cart2sph(N_recovery(1), N_recovery(2), N_recovery(3));
% rec_pt = [rec_az rec_el];
% for i=1:size(sh_c2,2)
%     sh_c_recov(i) = inverseSHT(sh_c2(:,i), rec_pt, 'real');
% end
% 
% % compare it with actual computed sh from a normal
% sh_recov_actual = find_blinnp_sh_coeff(N_recovery,w_o,tdesign_grid,n,Nord);
% 
% % first coeff test
% %recovered
% figure;imagesc(sh2sph(sh_c2(:,2),[100 100]))
% 
% %calculated
% sz = [100 100];
% [lon, lat] = meshgrid(linspace(-pi,pi,sz(2))', linspace(0,pi,sz(1))');
% sampling_grid = aziElev2aziIncl([lon(:) lat(:)]);
% temp = zeros(size(sampling_grid,1),1);
% for i=1:size(sampling_grid,1)
%     [Nix,Niy,Niz] = sph2cart(sampling_grid(i,1),sampling_grid(i,2),1);
%     Ni = [Nix;Niy;Niz];
%     t = find_blinnp_sh_coeff(Ni,w_o,tdesign_grid,n,1);
%     temp(i) = t(2);
% end
% temp = reshape(temp,sz(1),sz(2));
% figure;imagesc(temp)  
%     
%     
%     
% % real image test
% [rec_pts_az,rec_pts_el] = cart2sph(N_v(1,:),N_v(2,:),N_v(3,:));
% rec_pts = aziElev2aziIncl([rec_pts_az' rec_pts_el']);
% for i=1:size(sh_c2,2)
%     sh_c_recov_all(i,:) = inverseSHT(sh_c2(:,i), rec_pts, 'real');
% end
% 
% E_lm = repmat(const,1,size(sh_c_recov_all,2)).*sh_c_recov_all;
% Light_sh(121)=0;
% Y = getSH(Nord, rec_pts, 'real');
% 
% im = dot(  E_lm.*repmat(L_sh2',1,size(E_lm,2)),Y'   );
% im = reshape(im,size(N,1),size(N,2));
% figure;imshow(D + repmat(im*60/pi/2,1,1,3))
% 



%% SH specular - approach 2
%
R_v = 2*repmat(dot(repmat(w_o,1,size(N_v,2)),N_v),3,1).*N_v-repmat(w_o,1,size(N_v,2));
[rot_pts_az,rot_pts_el] = cart2sph(R_v(1,:),R_v(2,:),R_v(3,:));
rot_pts = aziElev2aziIncl([rot_pts_az' rot_pts_el']);

Nord = 10;
l = zeros((Nord+1)^2,1);
for i=0:Nord
    l(i^2+1:(i+1)^2) = i;
end
l = repmat(l,1,numel(n));
n = repmat(n(:)',size(l,1),1);
A_lm = exp(-l.^2/2./n);

E_lm = repmat(Light_sh,1,size(A_lm,2)).*A_lm;
S_v = rho_s'.*inverseSHT_my(E_lm, rot_pts, 'real');
S = zeros(size(is_face));
S(is_face) = S_v;
% S reshros(size(is_face));
% Seis_face)= S_s.*10;
% 
% [~, tdesign_grid] = getTdesign(2*Nord);
% tdesign_grid = aziElev2aziIncl(tdesign_grid); % convert to azi-incl
% 
% aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];
% 
% 
% A = max(cos(tdesign_grid(:,2)),0);
% A_lm = directSHT(Nord, A, tdesign_grid, 'real', []);
% 
% l = zeros((Nord+1)^2,1);
% for i=0:Nord
%     l(i^2+1:(i+1)^2) = i;
% end
% const = sqrt(4*pi./(2*l+1));
% 
% E_lm = const.*A_lm;
% for i=0:Nord
%     E_lm(i^2+1:(i+1)^2) = E_lm(median([i^2+1 (i+1)^2]));
% end
% 
% grid_sph = aziElev2aziIncl(tdesign_grid);
% % generate sh coeff for second term in the integral of blinn-phong. Do so
% % at each point of our grid
% Nord = 10;
% 
% % w_o(3) = w_o(3);
% % Light_sh(2:4) = Light_sh([3 4 2]);
% % Light_sh(121)=0;
% % probe env map at required points
% tdesign_grid_l = grid_sph;
% % tdesign_grid_l(:,2) = -tdesign_grid_l(:,2);
% tdesign_grid_l = aziElev2aziIncl(tdesign_grid_l);
% 
% Light_wi = inverseSHT(Light_sh, tdesign_grid_l, 'real');
% [wix_all,wiy_all,wiz_all] = sph2cart(grid_sph(:,1),grid_sph(:,2),1);
% wi_all = [wix_all';wiy_all';wiz_all'];
% H_all = (repmat(w_o ,1,size(wi_all,2))+wi_all);
% H_all = H_all./repmat(sqrt(sum(H_all.^2)),3,1);
% for i=1:size(grid_sph,1)
%     [Nix,Niy,Niz] = sph2cart(grid_sph(i,1),grid_sph(i,2),1);
%     Ni = [Nix;Niy;Niz];
%     N_dot_H = dot(repmat(Ni,1,size(H_all,2)),H_all);
%     %N_dot_wi = dot(repmat(Ni,1,size(wi_all,2)),wi_all);
%     A = rho_s*(n+2)/(2*pi)* Light_wi'.*(N_dot_H.^n);
%     %A = Light_wi';
%     A_sh = directSHT(Nord, A', tdesign_grid, 'real', []);
%     L(i,1) = inverseSHT(A_sh.*E_lm,tdesign_grid(i,:),'real');
% end
% L_sh = directSHT(Nord, L, tdesign_grid, 'real', []);
% 
% S = inverseSHT(L_sh, rec_pts, 'real');
% (S,size(N,1),size(N,2));
% figure;imshow(max(D,0) + max(repmat(S,1,1,3),0))

%% Subsurface diffuse
% i_o = [0 1  0 -1 0];
% j_o = [1 0 -1  0 0];
% % i_o = 0;
% % j_o = 0;
% 
% for clr=1:3
%     bssdf_r = bssdf(N_wi_theta,N_wo_theta,r,clr);
%     bssdf_r2r = bssdf(N_wi_theta,N_wo_theta,sqrt(2)*r,clr);
%     filtr = [0 1 0;1 0 1;0 1 0];
%     filtr2r = [1 0 1;0 0 0;1 0 1];
%     D(:,:,clr) = (conv2(bssdf_r,filtr,'same')+conv2(bssdf_r2r,filtr2r,'same')+bssdf(N_wi_theta,N_wo_theta,0,clr))*9*r^2;
% end
%% Naive diffuse
% naive_D = N_dot_wi*L_mag;
% naive_D(isnan(N(:,:,1))) = 0;
%% Output
% figure;imshow(repmat(naive_D,1,1,3).*alb)
% figure;imshow(D.*alb);
% 
% figure;imshow(repmat(S,1,1,3) + repmat(naive_D,1,1,3)/1.4)
% figure;imshow(repmat(S,1,1,3) + D);
% figure;imshow(max(D,0).*alb + max(repmat(S,1,1,3),0))
I = repmat(S,1,1,3) + alb.*D;

I_D = D;
end

