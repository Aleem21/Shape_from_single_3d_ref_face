function [ transmitance,theta_t ] = fresnel_trans( theta,n1,n2,type )
%FRESNEL_TRANS Summary of this function goes here
%   Detailed explanation goes here
if type==1
    theta_i = theta;
    theta_t = asin(n1/n2*sin(theta_i));
else
    theta_t = theta;
    theta_i = asin(sin(theta_t)*n2/n1);
end

t_parr = 2*sin(theta_t).*cos(theta_i)./(sin(theta_i+theta_t).*cos(theta_i-theta_t));
t_perp = 2*sin(theta_t).*cos(theta_i)./sin(theta_i+theta_t);
t = (t_parr+t_perp)/2;
transmitance = t.^2 .* n2.*cos(theta_t)/n1./cos(theta_i);

end

