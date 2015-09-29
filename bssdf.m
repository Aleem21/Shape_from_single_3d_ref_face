function [ f ] = bssdf( theta_i,theta_o,r,c )
%BSSDF Summary of this function goes here
%   Detailed explanation goes here
n_air = 1;
n_skin = 1.38;

if c == 1
    delta_a = 0.05;
    delta_s1 = 1.7;
elseif c==2
    delta_a = 0.24;
    delta_s1 = 1.36;
else
    delta_a = 0.43;
    delta_s1 = 1.7;
end

f = fresnel_trans(theta_i,n_air,n_skin,1)*...
    subsurf_profile(r,delta_a,delta_s1,n_skin)*...
    fresnel_trans(theta_o,n_skin,n_air,2);
end

