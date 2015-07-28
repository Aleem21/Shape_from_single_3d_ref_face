function [ cost, jacobian] = depth_cost_nonlin( z,i_p,i_q,i_bx,i_by,ncx,ncy,iz_reg,im,rhs_reg,l,rho,gaussVec)
%DEPTH_COST_NONLIL Summary of this function goes here
%   Detailed explanation goes here
%% data cost
p = z(i_p(:,1))-z(i_p(:,2));
q = z(i_q(:,1))-z(i_q(:,2));
cost_data = rho*l(1) + rho./(p.^2+q.^2+1).^0.5 .* (l(2)*p + l(3)*q - l(4))-im;

%% boundary conditions
cost_bound = (z(i_bx(:,1))-z(i_bx(:,2))).*ncx + (z(i_by(:,1))-z(i_by(:,2))).*ncy;

%% regularization cost
cost_reg = sum(z(iz_reg).*repmat(gaussVec',size(iz_reg,1),1),2)-rhs_reg;

%% sum it all up
% cost = sum(sum(cost_data.^2) + sum(cost_bound.^2) + sum(cost_reg.^2));
cost = [cost_data; cost_bound; cost_reg];


%% jacobian
if nargout >1
    nR = numel(cost);
    nC = numel(z);
    jacobian = sparse([],[],[],nR,nC,...
        numel(cost_data)*3+numel(cost_bound)*4+numel(cost_reg)*9);
    % data term
    constNumber = repmat(1:size(i_p,1),3,1)';
    d2 = (p.^2+q.^2+1);
    n = l(2)*p + l(3)*q - l(4);
    jacobian(sub2ind([nR nC],constNumber,[i_p i_q(:,1)])) =...
        repmat(rho./(d2.^1.5),1,3).* [l(2).*d2-n.*p -(l(2)+l(3)).*d2+n.*(p+q) l(3).*d2-n.*q];
    % boundary condition
    offset = size(i_p,1);
    constNumber = repmat(1:size(i_bx,1),1,1)' + offset;
    jacobian(sub2ind([nR nC],constNumber,i_bx(:,1))) = ncx;
    jacobian(sub2ind([nR nC],constNumber,i_by(:,1))) = ...
        jacobian(sub2ind([nR nC],constNumber,i_by(:,1))) + ncy;
    
    jacobian(sub2ind([nR nC],constNumber,i_bx(:,2))) = ...
        jacobian(sub2ind([nR nC],constNumber,i_bx(:,2))) - ncx;
    jacobian(sub2ind([nR nC],constNumber,i_by(:,2))) = ...
        jacobian(sub2ind([nR nC],constNumber,i_by(:,2))) - ncy;
    % reg term
    offset = offset + size(i_bx,1);
    constNumber = repmat(1:size(iz_reg,1),9,1)' + offset;
    jacobian(sub2ind([nR nC],constNumber,iz_reg)) = repmat(gaussVec',size(iz_reg,1),1);
    
end
end

