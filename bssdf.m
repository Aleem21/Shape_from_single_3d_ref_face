function [ f,f_filt ] = bssdf( theta_i,theta_o,r,c,delta_a,delta_s1 )
%BSSDF Summary of this function goes here
%   Detailed explanation goes here
import SHT.*;
n_air = 1;
n_skin = 1.38;
if nargin<6
    if c == 1
        delta_a = 0.12;
        delta_s1 = 0.7;
    elseif c==2
        delta_a = 0.34;
        delta_s1 = 0.6;
    else
        delta_a = 0.63;
        delta_s1 = 0.65;
    end
else
    delta_a = delta_a(c);
    delta_s1 = delta_s1(c);
end
if isempty(theta_i)
    f = [];
else
    f = fresnel_trans(theta_i,n_air,n_skin,1).*...
        subsurf_profile(r,delta_a,delta_s1,n_skin).*...
        fresnel_trans(theta_o,n_skin,n_air,2)/pi;
end
if nargout>1
    dist = zeros(5);
    mid = [3,3];
    [y,x] = meshgrid(1:5,1:5);
    xmid = dist+mid(1);
    ymid = dist+mid(2);
    dist = sqrt((x-xmid).^2 + (y-ymid).^2);
    f_filt=subsurf_profile(dist*r,delta_a,delta_s1,n_skin);
end
end

