function [ n_mean,n_var,rho_mean,rho_var ] = load_merl_shading( merl_path )
%LOAD_MERL_ALB Summary of this function goes here
%   Detailed explanation goes here
load([merl_path '\bp_data\bp_data.mat'])
[n_mean,n_var]=compute_gaussian_mixture(bp(:,1),bp_region);
[rho_mean,rho_var]=compute_gaussian_mixture(bp(:,2)*50,bp_region);

end

