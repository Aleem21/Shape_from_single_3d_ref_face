function [ rendered,rendered_c,z,scales ] = read_render_USF( impath,Rpose,sz,is_full)
%READ_RENDER_USF Summary of this function goes here
%   Detailed explanation goes here
if nargin <3
    sz = [300 300];
end
if nargin<4
    is_full = 0;
end
[pts,tri,rgb] = read_USF_eko(impath,512,512,0,0,is_full);
pts_rotated = Rpose*[pts; ones(1,size(pts,2))];
pts_rotated = pts_rotated(1:3,:);


[rendered_c,z] = render_rgb_USF(pts_rotated,tri,rgb,sz(1),sz(2));
rendered_c = double(im2double(rendered_c));
rendered = rgb2gray(rendered_c);
z = double(z)*min(sz)/2;
scales = [1/9 1/9 1/9];
end