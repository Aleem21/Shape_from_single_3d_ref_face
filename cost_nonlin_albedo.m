function [ cost, jacobian] = cost_nonlin_albedo( albedo,p,q,iz_reg,im,rhs_reg,l,gaussVec,eye_mask,is_face)

%% data cost
dark_factor = eye_mask;
cost_data = (albedo*l(1) + albedo./(p.^2+q.^2+1).^0.5 .* (l(2)*p + l(3)*q - l(4))-im).*dark_factor;

%% regularization cost
cost_reg = sum(albedo(iz_reg).*repmat(gaussVec',size(iz_reg,1),1),2)-rhs_reg;

%% sum it all up
% cost = sum(sum(cost_data.^2) + sum(cost_bound.^2) + sum(cost_reg.^2));
cost = [cost_data; cost_reg];


%% jacobian
if nargout >1
    nR = numel(cost);
    nC = numel(albedo);
    d = (p.^2+q.^2+1).^0.5;
    n = l(2)*p + l(3)*q - l(4);
    % data term
    constNumber1 = 1:numel(albedo);
    data_rhs = (l(1)+n./d).*dark_factor;
    
    % reg term
    offset = numel(cost_data);
    constNumber2 = repmat(1:size(iz_reg,1),9,1)' + offset;
    reg_rhs = repmat(gaussVec',size(iz_reg,1),1);
    jacobian = sparse([constNumber1(:); constNumber2(:)]...
        ,[constNumber1(:); iz_reg(:)]...
        ,[data_rhs(:); reg_rhs(:)],...
        nR,nC);
end

end


