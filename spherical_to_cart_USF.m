function [ x,y,z ] = spherical_to_cart_USF( spherical,r,c )
%SPHERICAL_TO_CART_USF convert spherical coordinates output in USF to
%cartesian
spherical = spherical /2^15 * 524272e-6;
theta = repmat(linspace(0,2*pi,c),r,1);
x = spherical.* cos(theta);
y = spherical.* sin(theta);
z = repmat(linspace(-615e-6*512/2,615e-6*512/2,r)',1,c);
end

