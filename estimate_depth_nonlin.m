function [ depth,alb_out,in_face,l_sh] = estimate_depth_nonlin...
    ( alb_ref, im, z_ref, sh_coeff,lambda1,lambda2,lambda_bound,max_iter,bound_type,...
    jack,eye_mask,is_dz_depth,is_l_sh,lambda_dz,lambda_reg2,z_gnd,algo,talk)
%ESTIMATE_DEPTH Summary of this function goes here
%   Detailed explanation goes here
if nargin <15
    talk = 1;
end
if nargin <14
    algo = 'levenberg-marquardt';
end
if nargin<13
    is_gnd = 0;
else
    is_gnd = 1;
    if isempty(z_gnd)
        is_gnd = 0;
    end
end
%% Prepocessing
[ face,face_inds, inface_inds,in_face] = preprocess_estimate_depth( z_ref );

%% Optimization - depth
[ costfun_z, is_face,nData,nBound,nReg]...
    = get_depth_costfun( z_ref, im, alb_ref, sh_coeff,eye_mask, lambda1,...
    lambda_bound,bound_type,is_dz_depth,is_l_sh,lambda_dz,lambda_reg2);

init_z = [z_ref(is_face); sh_coeff(1:4)];
% options = optimset('Display','iter-detailed','maxIter',100,'JacobPattern',jacobianPattern);
% options = optimset('Display','iter-detailed','maxIter',200,...
%     'Jacobian','on','JacobMult',@jacobMultFnc,'JacobPattern',jacobianPattern);
options = optimset('Display','iter-detailed','maxIter',max_iter,...
    'Jacobian',jack,'Algorithm',algo,'maxFunEvals',Inf,'TolX',0.5); 

% options = optimset('Display','iter-detailed','maxIter',max_iter,...
%     'JacobPattern',jacobianPattern); %,'Algorithm','levenberg-marquardt'
% options = optimset('maxIter',1,'DerivativeCheck','on','Jacobian','on');
% [z]=lsqnonlin(costfun_z,init_z,[],[],options);

z=lsqnonlin(costfun_z,init_z,[],[],options);
% options = optimset('Display','iter-detailed','maxIter',max_iter,...
%     'Jacobian','off','Algorithm',algo,'maxFunEvals',Inf); %
% z2=lsqnonlin(costfun_z,z,[],[],options);
l_sh = z(end-3:end);
z = z(1:end-4);
%% Post precessing
depth = NaN(size(alb_ref));
depth(face_inds) = z;
depth(~face) = NaN;

offset = mean(depth(inface_inds)) - mean(z_ref(inface_inds));
depth = depth-offset;

%% Optimization - albedo
if nargout>1
    init_alb = alb_ref(is_face);
    [ costfun_alb,jacobianPattern_alb ] = get_albedo_costfun( depth, im,alb_ref, sh_coeff, eye_mask,lambda2);
    options = optimset('Display','iter-detailed','maxIter',max_iter,...
        'Jacobian',jack,'JacobPattern',jacobianPattern_alb,'Algorithm','levenberg-marquardt'); %
    alb_est = lsqnonlin(costfun_alb,init_alb,[],[],options);
    alb_out = zeros(size(alb_ref));
    alb_out(is_face) = alb_est;
    alb_out(~in_face)=0;
end
depth(~in_face)=nan;
%% showing output
if talk
    cost_ref = costfun_z(init_z).^2;
    cost_est = costfun_z([z(:); l_sh(:)]).^2;
    if is_gnd
        z_gndd = z_gnd(is_face);
        cost_gt= costfun_z(z_gndd).^2;
    end
    
    s_data = 1;
    e_data = nData;
    s_bnd = e_data+1;
    e_bnd = s_bnd-1 + nBound;
    s_reg = e_bnd+1;
    e_reg = s_reg-1 + nReg;
    s_reg2 = e_reg+1;
    e_reg2 = numel(cost_ref);
    
    figure;
    subplot(2,1,1);
    plot(s_data:e_data,cost_ref(s_data:e_data),...
        s_bnd:e_bnd,cost_ref(s_bnd:e_bnd),...
        s_reg:e_reg,cost_ref(s_reg:e_reg),...
        s_reg2:e_reg2,cost_ref(s_reg2:e_reg2));
    title('Reference cost')
    subplot(2,1,2);
    plot(s_data:e_data,cost_est(s_data:e_data),...
        s_bnd:e_bnd,cost_est(s_bnd:e_bnd),...
        s_reg:e_reg,cost_est(s_reg:e_reg),...
        s_reg2:e_reg2,cost_est(s_reg2:e_reg2));
    title('Optimized cost');
    
    if is_gnd
        figure;
        plot(s_data:e_data,cost_gt(s_data:e_data),...
            s_bnd:e_bnd,cost_gt(s_bnd:e_bnd),...
            s_reg:e_reg,cost_gt(s_reg:e_reg));
        title('Ground truth cost')
    end

end
