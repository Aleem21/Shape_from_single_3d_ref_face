function [ x,y,z ] = spherical_to_cart_USF( spherical,r,c )
%SPHERICAL_TO_CART_USF convert spherical coordinates output in USF to
%cartesian

theta = repmat(linspace(0,2*pi,c),r,1);
x = spherical.* cos(theta);
y = spherical.* sin(theta);
z = repmat(linspace(-9000,9000,r)',1,c);
end

