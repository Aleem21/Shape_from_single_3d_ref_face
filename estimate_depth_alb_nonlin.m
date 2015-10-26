function [ depth,alb_out,is_face ] = estimate_depth_alb_nonlin...
    ( alb_ref, im, z_ref, sh_coeff,lambda1,lambda2,lambda_bound,max_iter,...
    bound_type,jack,eye_mask,is_alb_dz,z_gnd,algo,talk)
%ESTIMATE_DEPTH Summary of this function goes here
%   Detailed explanation goes here
if nargin <15
    talk = 1;
end
if nargin < 14
    algo = 'levenberg-marquardt';
end
if nargin<13
    is_gnd = 0;
else
    is_gnd = 1;
end
%% Prepocessing

[ face,face_inds, inface_inds] = preprocess_estimate_depth( z_ref );

%% Optimization
[ costfun_z_alb, is_face,nData,nBound,nReg,jacobianPattern_z ]...
    = get_depth_alb_costfun( z_ref, im, alb_ref, sh_coeff,eye_mask, lambda1,lambda2,lambda_bound,bound_type,is_alb_dz);
init_z = z_ref(is_face);
init_alb = alb_ref(is_face);
init_z_alb = [init_z ; init_alb];
% options = optimset('Display','iter-detailed','maxIter',100,'JacobPattern',jacobianPattern);
% options = optimset('Display','iter-detailed','maxIter',200,...
%     'Jacobian','on','JacobMult',@jacobMultFnc,'JacobPattern',jacobianPattern);

options = optimset('Display','iter-detailed','maxIter',max_iter,...
    'Jacobian',jack,'JacobPattern',jacobianPattern_z,'Algorithm',algo); %
% options = optimset('Display','iter-detailed','maxIter',max_iter,...
%     'JacobPattern',jacobianPattern); %,'Algorithm','levenberg-marquardt'
% options = optimset('maxIter',1,'DerivativeCheck','on','Jacobian','on');
% [z]=lsqnonlin(costfun,init_z,[],[],options);

z_alb=lsqnonlin(costfun_z_alb,init_z_alb,[],[],options);

%% Post precessing
z = z_alb(1:end/2);
alb_est = z_alb(end/2+1:end);
alb_out = zeros(size(alb_ref));
alb_out(is_face) = alb_est;

depth = NaN(size(alb_ref));
depth(face_inds) = z;
depth(~face) = NaN;

offset = mean(depth(inface_inds)) - mean(z_ref(inface_inds));
depth = depth-offset;
cost_ref = costfun_z_alb(init_z_alb).^2;
cost_est = costfun_z_alb(z_alb).^2;
figure;plot(cost_ref);hold on;
plot(cost_est);
%% showing output
if talk
    %     cost_ref = costfun_z(init_z).^2;
    %     cost_est = costfun_z(z).^2;
    %     if is_gnd
    %         z_gndd = z_gnd(is_face);
    %         cost_gt= costfun_z(z_gndd).^2;
    %     end
    %
    %     s_data = 1;
    %     e_data = nData;
    %     s_bnd = e_data+1;
    %     e_bnd = s_bnd-1 + nBound;
    %     s_reg = e_bnd+1;
    %     e_reg = numel(cost_ref);
    %
    %     figure;
    %     subplot(2,1,1);
    %     plot(s_data:e_data,cost_ref(s_data:e_data),...
    %         s_bnd:e_bnd,cost_ref(s_bnd:e_bnd),...
    %         s_reg:e_reg,cost_ref(s_reg:e_reg));
    %     title('Reference cost')
    %     subplot(2,1,2);
    %     plot(s_data:e_data,cost_est(s_data:e_data),...
    %         s_bnd:e_bnd,cost_est(s_bnd:e_bnd),...
    %         s_reg:e_reg,cost_est(s_reg:e_reg));
    %     title('Optimized cost');
    %
    %     if is_gnd
    %         figure;
    %         plot(s_data:e_data,cost_gt(s_data:e_data),...
    %             s_bnd:e_bnd,cost_gt(s_bnd:e_bnd),...
    %             s_reg:e_reg,cost_gt(s_reg:e_reg));
    %         title('Ground truth cost')
    %     end
end
