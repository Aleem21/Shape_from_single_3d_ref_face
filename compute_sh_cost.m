function [ error ] = compute_sh_cost( im, n, alb,l, nComp )
%COMPUTE_SH_COST Summary of this function goes here
%   Detailed explanation goes here
a0 = pi;
a1 = 2*pi/sqrt(3);
a2= 2*pi/sqrt(8);

c0 =    sqrt(1  /(4*pi));
c1 =    sqrt(3  /(4*pi));
c2 = 3* sqrt(5  /(12*pi));

nx = n(:,:,1);
nx = -nx(:)';        %change to row vector
ny = n(:,:,2);
ny = -ny(:)';
nz = n(:,:,3);
nz = nz(:)';

Y = [   ones(size(nx));
        nx;
        ny;
        nz;
        nx.*ny;
        nx.*nz;
        ny.*nz;
        (nx.^2 - ny.^2);
        (3*nz.^2-1) ];

% pick 1st 4 coefficients and remove nan values

Y = Y(1:nComp,~isnan(nx));
%vectorize image as row vector and remove nan values
im = im(:);
im = im(~isnan(nx))';
alb = alb(:);
alb = alb(~isnan(nx))';

bad = im<0.05 | im>0.95 | alb< 0.05;

im(bad) = [];
alb(bad) = [];
Y(:,bad) = [];
%account for brdf
rho = ones(size(im));
radiosity = im./rho./alb;
l = l.* [c0; c1; c1; c1;  c2; c2; c2; c2/2; c2/(12)^.5 ];

error = sum((radiosity - (l(1:nComp)'*Y)).^2);


end

