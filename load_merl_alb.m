function [ albedo_mean,var ] = load_merl_alb( merl_path )
%LOAD_MERL_ALB Summary of this function goes here
%   Detailed explanation goes here
load([merl_path 'ts_data.mat'])
[albedo_mean,var]=compute_gaussian_mixture(albedo,albedo_region);

end

