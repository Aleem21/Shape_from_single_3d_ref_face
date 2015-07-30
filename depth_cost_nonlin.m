function [ cost, jacobian] = depth_cost_nonlin( z,i_p,i_q,i_bx,i_by,ncx,ncy,iz_reg,im,rhs_reg,l,rho,gaussVec,type,i_bound,val_bound)
%DEPTH_COST_NONLIL Summary of this function goes here
%   Detailed explanation goes here
%% data cost
p = z(i_p(:,1))-z(i_p(:,2));
q = z(i_q(:,1))-z(i_q(:,2));
cost_data = rho*l(1) + rho./(p.^2+q.^2+1).^0.5 .* (l(2)*p + l(3)*q - l(4))-im;

%% boundary conditions
if type==1
    cost_bound = (z(i_bx(:,1))-z(i_bx(:,2))).*ncx + (z(i_by(:,1))-z(i_by(:,2))).*ncy;
else
    cost_bound = sum(z(i_bound).*val_bound,2);
end
%% regularization cost
cost_reg = sum(z(iz_reg).*repmat(gaussVec',size(iz_reg,1),1),2)-rhs_reg;

%% sum it all up
% cost = sum(sum(cost_data.^2) + sum(cost_bound.^2) + sum(cost_reg.^2));
cost = [cost_data; cost_bound; cost_reg];


%% jacobian
if nargout >1
    nR = numel(cost);
    nC = numel(z);
    d2 = (p.^2+q.^2+1);
    n = l(2)*p + l(3)*q - l(4);
    % data term
    constNumber1 = repmat(1:size(i_p,1),3,1)';
    data_rhs = repmat(rho./(d2.^1.5),1,3).*...
        [l(2).*d2-n.*p -(l(2)+l(3)).*d2+n.*(p+q) l(3).*d2-n.*q];
    
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
    if type==1
        jacobian = sparse([constNumber1(:); constNumber2(:); constNumber3(:)]...
        ,[i_p(:); i_q(:,1);i_bx(:);i_by(:);iz_reg(:)]...
        ,[data_rhs(:); bnd_rhs(:); reg_rhs(:)],...
        nR,nC);
    else
    jacobian = sparse([constNumber1(:); constNumber2(:); constNumber3(:)]...
        ,[i_p(:); i_q(:,1);i_bnd;i_by(:);iz_reg(:); ]...
        ,[data_rhs(:); bnd_rhs(:); reg_rhs(:)],...
        nR,nC);
    end
end
end

