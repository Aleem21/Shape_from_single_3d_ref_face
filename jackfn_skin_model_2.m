function [ cost,J ] = jackfn_skin_model_2(params,costfn,costfn_2,options)
%COSTFN_SKIN_MODEL Summary of this function goes here
%   Detailed explanation goes here
cost = costfn(params);
if nargout >1
    
    [J,fevals] = compute_Jacobian(params,costfn_2,options);
end

% cost = cost*0;
% cost(27495)=1-params(13748);
end

