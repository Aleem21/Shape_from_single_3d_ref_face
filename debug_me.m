function [ d,lhs,rhs ] = debug_me( ind,indr,indb,z,r_face,c_face,im,alb,N,l0,l1,l2,l3 )
%DEBUG_ME Summary of this function goes here
%   Detailed explanation goes here
r = r_face([ind indr indb]);
c = c_face([ind indr indb]);

q = z(r(2),c(2))-z(r(1),c(1));
p = z(r(3),c(3))-z(r(1),c(1));
I = im(r(1),c(1));
lhs = I - alb(r(1),c(1))*l0 + alb(r(1),c(1))/N(r(1),c(1))*l3;
rhs = alb(r(1),c(1))/N(r(1),c(1))*(l1*p+l2*q);
d = lhs-rhs;
fprintf('r %d, c %d',r(1),c(1));
end

