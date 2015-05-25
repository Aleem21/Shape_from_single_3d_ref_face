function [ im ] = render_model_noGL( n, sh_coeff, alb,talk )
%RENDER_MODEL Summary of this function goes here
%   Detailed explanation goes here
if nargin<4
    talk = 0;
end
sh_coeff(10) = 0;
sh_coeff(10) = [];
alb = alb(:)';
nx = n(:,:,1);
nx = nx(:)';        %change to row vector
ny = n(:,:,2);
ny = ny(:)';
nz = n(:,:,3);
nz = nz(:)';

Y = [ones(size(nx));
    nx;
    ny;
    nz;
    nx.*ny;
    nx.*nz;
    ny.*nz;
    (nx.^2 - ny.^2);
    (3*nz.^2-1) ];
im = alb.*(sh_coeff'*Y);
im(im<0) = 0;
im = reshape(im,size(n,1),size(n,2));

if talk
    figure; imshow(im);
end
end

