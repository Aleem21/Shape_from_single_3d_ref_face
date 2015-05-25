function [ l ] = estimate_lighting( n, alb, im,nComp,talk )
%ESTIMATE_LIGHTING etimate lighting spherical harmonic coefficients given
%the normal map and input image
%
%   im must be grayscale
%   n must be of dimensions same as image with 3rd dimension having 3
%   slices containing nx,ny,nz

if nargin<4
    nComp = 4;
end
if nargin<5
    talk = 0;
end


if talk
    rad = im./alb;
    rad(abs(rad)>1.5)=nan;
    figure;imagesc(rad);
end

a0 = pi;
a1 = 2*pi/sqrt(3);
a2= 2*pi/sqrt(8);

c0 =    sqrt(1  /(4*pi));
c1 =    sqrt(3  /(4*pi));
c2 = 3* sqrt(5  /(12*pi));

nx = n(:,:,1);
nx = nx(:)';        %change to row vector
ny = n(:,:,2);
ny = ny(:)';
nz = n(:,:,3);
nz = nz(:)';

Y = [   c0*     ones(size(nx));
        c1*     nx;
        c1*     ny;
        c1*     nz;
        c2*     nx.*ny;
        c2*     nx.*nz;
        c2*     ny.*nz;
       c2/2*    (nx.^2 - ny.^2);
  c2/(12)^.5*   (3*nz.^2-1) ];

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

l = radiosity * pinv(Y);
l(10) = 0;
l(10) = [];
l = l'.* [c0; c1; c1; c1;  c2; c2; c2; c2/2; c2/(12)^.5 ];
end

