function [ rendered,rendered_c,z ] = read_render_USF( impath,Rpose )
%READ_RENDER_USF Summary of this function goes here
%   Detailed explanation goes here
[pts,tri,rgb] = read_USF_eko(impath,512,512);
pts_rotated = Rpose*[pts; ones(1,size(pts,2))];
pts_rotated = pts_rotated(1:3,:);


[rendered_c,z] = render_rgb_USF(pts_rotated,tri,rgb,300,300);
rendered_c = double(im2double(rendered_c));
rendered = rgb2gray(rendered_c);
z = double(z);
end

