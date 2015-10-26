function [ Pd] = subsurf_profile( r,delta_a,delta_s1,n )
%SUBSURF_PROFILE Summary of this function goes here
%   Detailed explanation goes here


F_dr = -1.440/n^2 + 0.710/n + 0.668 + 0.0636*n;% The Magic of Computer 
                                               % Graphics By Noriko Kurachi
                                               
J = (1+F_dr)/(1-F_dr);

delta_t1 = delta_a+delta_s1;

z_r = 1/delta_t1;
z_v = (1+4*J/3)/delta_t1;

d_r = sqrt(r.^2+z_r^2);
d_v = sqrt(r.^2+z_v^2);

alpha_1 = delta_s1/delta_t1;

delta_tr = sqrt(3*delta_a*delta_t1);
Pd = alpha_1*z_r * (1+delta_tr*d_r).*exp(-delta_tr*d_r)./(4*pi*d_r.^3) + ...
     alpha_1*z_v * (1+delta_tr*d_v).*exp(-delta_tr*d_v)./(4*pi*d_v.^3);

end

