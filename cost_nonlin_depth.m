function [ cost, jacobian] = cost_nonlin_depth( z_l,i_p,i_q,i_bx,i_by,ncx,ncy,iz_reg,im,rhs_reg,lambda_dz,lambda_reg2,l,rho,z_ref,gaussVec,type,is_dz_depth,is_l_sh,eye_mask,i_bound,val_bound,is_face,face)
%DEPTH_COST_NONLIL Summary of this function goes here
%   Detailed explanation goes here
%% data cost
z = z_l(1:end-4);
l = z_l(end-3:end);
p = z(i_p(:,1))-z(i_p(:,2));
q = z(i_q(:,1))-z(i_q(:,2));
dark_factor = eye_mask;
cost_data = (rho*l(1) + rho./(p.^2+q.^2+1).^0.5 .* (l(2)*p + l(3)*q - l(4))-im).*dark_factor;

synth = fill_var(rho*l(1) + rho./(p.^2+q.^2+1).^0.5 .* (l(2)*p + l(3)*q - l(4)),is_face,0);
im2 = fill_var(im,is_face,0);
z_synth = fill_var(z_l(i_p(:,1)),is_face,nan);
% p_synth = fill_var(p,is_face,0);
% q_synth = fill_var(q,is_face,0);
 


%% boundary conditions
if type==1
    cost_bound = (z(i_bx(:,1))-z(i_bx(:,2))).*ncx + (z(i_by(:,1))-z(i_by(:,2))).*ncy;
else
    cost_bound = sum(z(i_bound).*val_bound,2);
end
%% regularization cost
cost_reg = sum(z(iz_reg).*repmat(gaussVec',size(iz_reg,1),1),2)-rhs_reg;
cost_reg_im = fill_var(cost_reg,is_face,0);
if is_dz_depth
    cost_reg_dz = (z(iz_reg(:,5))-z_ref)*lambda_dz;
else
    cost_reg_dz = [];
end
lap_op = [0.25 0.5 0.25;0.5 -3 0.5;0.25 0.5 0.25]*lambda_reg2;
cost_reg2 = sum(  z(iz_reg).*repmat(lap_op(:)',size(iz_reg,1),1)  ,2  );
%% sum it all up
% cost = sum(sum(cost_data.^2) + sum(cost_bound.^2) + sum(cost_reg.^2));
cost = [cost_data; cost_bound; cost_reg; cost_reg_dz; cost_reg2];

reg_synth = fill_var(cost_reg,is_face,0);
data_synth = fill_var(cost_data,is_face,0);

%% jacobian
if nargout >1
    nR = numel(cost);
    nC = numel(z)+4;
    d2 = (p.^2+q.^2+1);
    n = l(2)*p + l(3)*q - l(4);
    % data term
    constNumber1 = repmat(1:size(i_p,1),4,1)';
    data_rhs = repmat(rho./(d2.^1.5),1,4).*...
        [l(2).*d2-n.*p  -l(2).*d2+n.*p  l(3).*d2-n.*q -l(3).*d2+n.*q].*repmat(dark_factor,1,4);
    
    % boundary condition
    
    
    
    offset = size(i_p,1);
    if type==1
        constNumber2 = repmat(1:size(i_bx,1),4,1)' + offset;
        bnd_rhs = [ncx; -ncx; ncy; -ncy];
    else
        bnd_rhs = [];
        i_bnd = [];
        constNumber2 = [];
        for i=1:size(val_bound,1);
            vals = val_bound(i,:);
            is = i_bound(i,:);
            [C,~,bi] = unique(is);
            cols = 1:numel(vals);
            rows = bi';
            temp = zeros(numel(C),numel(bi));
            temp(sub2ind(size(temp),rows,cols)) = vals;
            temp = sum(temp,2);
            i_bnd = [i_bnd ; C(:)];
            bnd_rhs = [bnd_rhs; temp(:)];
            constNumber2 = [constNumber2; ones(numel(temp),1)*i];
        end
        constNumber2 = constNumber2 + offset;      
    end
    % reg term
    offset2 = offset + numel(cost_bound);
    constNumber3 = repmat(1:size(iz_reg,1),9,1)' + offset2;
    reg_rhs = repmat(gaussVec',size(iz_reg,1),1);
    
    offset3 = offset2 + numel(cost_reg);
    if is_dz_depth
        constNumber4 = (1:size(iz_reg,1)) + offset3;
        col4 = iz_reg(:,5);
        rhs_reg_dz = ones(1,numel(col4))*lambda_dz;
    else
        constNumber4 = [];
        col4 = [];
        rhs_reg_dz = [];
    end
    offset4 = offset3 + numel(cost_reg_dz);
    constNumber5 = repmat(1:size(iz_reg,1),9,1)' + offset4;
    reg_rhs_2 = repmat(lap_op(:)',size(iz_reg,1),1);
    
    
    %spherical harmonic lighting
    constNumber6 = repmat(1:size(i_p,1),4,1);
    col6 = repmat((1:4)',1,size(i_p,1))+numel(z);
    sh_rhs = [rho [p q p*0-1].*repmat(rho./sqrt(d2),1,3)]';
    if ~is_l_sh
        sh_rhs = sh_rhs*0;
    end
    if type==1
        jacobian = sparse([constNumber1(:); constNumber2(:); constNumber3(:); constNumber4(:)]...
        ,[i_p(:); i_q(:);i_bx(:);i_by(:);iz_reg(:);col4(:)]...
        ,[data_rhs(:); bnd_rhs(:); reg_rhs(:); rhs_reg_dz(:)],...
        nR,nC);
    else
    jacobian = sparse([constNumber1(:); constNumber2(:); constNumber3(:);...
        constNumber4(:); constNumber5(:);constNumber6(:)]...
        ,[i_p(:); i_q(:);i_bnd;i_by(:);iz_reg(:);col4(:);iz_reg(:);col6(:)]...
        ,[data_rhs(:); bnd_rhs(:); reg_rhs(:); rhs_reg_dz(:); reg_rhs_2(:);sh_rhs(:)],...
        nR,nC);
    end
end
end
function out = fill_var(x,is_face,to_fill)
if nargin<3
    to_fill = 0;
end
out = zeros(size(is_face));
out(:) = to_fill;
out(is_face) = x;
end

