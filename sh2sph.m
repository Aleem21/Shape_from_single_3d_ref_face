function [image,x,y,z] = sh2sph(sh_coeff,sz,draw)
if nargin<3
    draw=0;
end

c0 =    sqrt(1  /(4*pi));
c1 =    sqrt(3  /(4*pi));
c2 = 3* sqrt(5  /(12*pi));

if size(sh_coeff,2)>1
    for i=1:size(sh_coeff,2)
        [image(:,:,i),x,y,z] = sh2sph(sh_coeff(:,i),sz,0);
    end
else
    sh_coeff(10) = 0;
    sh_coeff(10) = [];
    sh_coeff = sh_coeff ./ [c0 c1 c1 c1 c2 c2 c2 c2/2  c2/2/sqrt(3)]';
    sh_coeff([2 3 4 5 6 7 8 9]) = sh_coeff([3 4 2 5 7 9 6 8]);
    sh_coeff([4 5 8]) = -sh_coeff([4 5 8]);
    [lon, lat] = meshgrid(linspace(-pi,pi,sz(2))', linspace(pi/2,-pi/2,sz(1))');
    [x,y,z] = lla2ecef(lat,lon,ones(size(lat)));
    sampling_grid = [lon(:) pi/2-lat(:)];
    env_sh = inverseSHT(sh_coeff, sampling_grid, 'real');
    image = reshape(env_sh,sz);
    
end
if draw
    figure;surf(x,y,z,repmat(image,1,1,3),'edgealpha',0);
    axis equal;
    view(180,-90)
    colormap gray;
    xlabel('x');
    ylabel('y');
    zlabel('z');
end
end