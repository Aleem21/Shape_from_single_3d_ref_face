function [ I,S,I_D,naive_D ] = render_DS( N,L,w_o,n,rho_s,amb,im_c,r )
%RENDER_DS Summary of this function goes here
%   Detailed explanation goes here

%% Definitions
% w_o                       - viewing direction
% rho_s                     - specular scale
% n                         - specular exponent
%% Specular part
% light direction
L_d = L/norm(L);
L_mag = norm(L);
% half angle vector
H = -(L+w_o)/norm(L_d+w_o);
N_v(1,:) = reshape(N(:,:,1),1,[]);
N_v(2,:) = reshape(N(:,:,2),1,[]);
N_v(3,:) = reshape(N(:,:,3),1,[]);
H_v = repmat(H,1,size(N_v,2));
cos_delta = dot(N_v,H_v);
cos_delta = max(reshape(cos_delta,size(N,1),size(N,2)),0);

% S = zeros(size(N,1),size(N,2));
% for i=1:size(S,1)
%     for j=1:size(S,2)
%        S(i,j) = rho_s*(n+2)/(2*pi) * 
%     end
% end
N_dot_wi = dot(N_v,repmat(-L_d,1,size(N_v,2)));
N_dot_wi = max(reshape(N_dot_wi,size(N,1),size(N,2)),0);
S = rho_s * (n+2)/(2*pi) * cos_delta.^n * L_mag .* N_dot_wi;

%% Subsurface diffuse
i_o = [0 1  0 -1 0];
j_o = [1 0 -1  0 0];
% i_o = 0;
% j_o = 0;
N_wi_theta = acos(N_dot_wi);
N_wo_theta = acos(reshape(dot(N_v,repmat(-w_o,1,size(N_v,2)) ),size(N,1),size(N,2)));
D = zeros(size(N,1),size(N,2),1);
for clr=1:1
for i=2:size(N,1)-1
    disp(i)
    for j=2:size(N,2)-1
        for ind = 1:numel(i_o)
            D(i,j,clr) = D(i,j,clr) +...
                bssdf(N_wi_theta(i+i_o(ind),j+j_o(ind)),N_wo_theta(i+i_o(ind),j+j_o(ind)),r,1);
        end
    end
end
end
D = repmat(D,1,1,3);
%% Naive diffuse
naive_D = N_dot_wi*L_mag;
naive_D(isnan(N(:,:,1))) = 0;
%% Output
figure;imshow(repmat(naive_D,1,1,3)/1.4)
figure;imshow(2.5*D/1.2);

figure;imshow(repmat(S,1,1,3) + repmat(naive_D,1,1,3)/1.4)
figure;imshow(repmat(S,1,1,3) + 2.5*D/1.2);

I = repmat(S,1,1,3) + 2.5*D/1.2;

I_D = 2.5*D/1.2;

end

