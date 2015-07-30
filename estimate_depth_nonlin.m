function [ depth ] = estimate_depth_nonlin...
    ( alb_ref, im, z_ref, sh_coeff,lambda1,lambda_bound,max_iter,bound_type,jack,eye_mask,z_gnd,talk)
%ESTIMATE_DEPTH Summary of this function goes here
%   Detailed explanation goes here
if nargin <12
    talk = 1;
end
if nargin<11
    is_gnd = 0;
else
    is_gnd = 1;
end
%% Prepocessing

[ face,face_inds, inface_inds] = preprocess_estimate_depth( z_ref );

%% Optimization
[ costfun, is_face,nData,nBound,nReg,jacobianPattern ]...
    = get_costfun( z_ref, im, alb_ref, sh_coeff,eye_mask, lambda1,lambda_bound,bound_type);

init_z = z_ref(is_face);
% options = optimset('Display','iter-detailed','maxIter',100,'JacobPattern',jacobianPattern);
% options = optimset('Display','iter-detailed','maxIter',200,...
%     'Jacobian','on','JacobMult',@jacobMultFnc,'JacobPattern',jacobianPattern);
options = optimset('Display','iter-detailed','maxIter',max_iter,...
    'Jacobian',jack,'JacobPattern',jacobianPattern,'Algorithm','levenberg-marquardt'); %
% options = optimset('Display','iter-detailed','maxIter',max_iter,...
%     'JacobPattern',jacobianPattern); %,'Algorithm','levenberg-marquardt'
% options = optimset('maxIter',1,'DerivativeCheck','on','Jacobian','on');
% [z]=lsqnonlin(costfun,init_z,[],[],options);

[z,fval,residual,exitflag,output]=lsqnonlin(costfun,init_z,[],[],options);

%% Post precessing
depth = NaN(size(alb_ref));
depth(face_inds) = z;
depth(~face) = NaN;

offset = mean(depth(inface_inds)) - mean(z_ref(inface_inds));
depth = depth-offset;


%% showing output
if talk
    cost_ref = costfun(init_z).^2;
    cost_est = costfun(z).^2;
    if is_gnd
        z_gndd = z_gnd(is_face);
        cost_gt= costfun(z_gndd).^2;
    end
    
    s_data = 1;
    e_data = nData;
    s_bnd = e_data+1;
    e_bnd = s_bnd-1 + nBound;
    s_reg = e_bnd+1;
    e_reg = numel(cost_ref);
    
    figure;
    subplot(2,1,1);
    plot(s_data:e_data,cost_ref(s_data:e_data),...
        s_bnd:e_bnd,cost_ref(s_bnd:e_bnd),...
        s_reg:e_reg,cost_ref(s_reg:e_reg));
    title('Reference cost')
    subplot(2,1,2);
    plot(s_data:e_data,cost_est(s_data:e_data),...
        s_bnd:e_bnd,cost_est(s_bnd:e_bnd),...
        s_reg:e_reg,cost_est(s_reg:e_reg));
    title('Optimized cost');
    
    if is_gnd
        figure;
        plot(s_data:e_data,cost_gt(s_data:e_data),...
            s_bnd:e_bnd,cost_gt(s_bnd:e_bnd),...
            s_reg:e_reg,cost_gt(s_reg:e_reg));
        title('Ground truth cost')
    end
end
