function [ z,rendered ] = computer_depth_USF( pts,tri,rgb,xrange,yrange,im_c,talk )
%COMPUTER_DEPTH_USF Summary of this function goes here
%   Detailed explanation goes here
if nargin<4
    talk = 0;
end
% scale xrange space of model into [-1,1]
pts(1,:) = (pts(1,:)-xrange(1))/diff(xrange)*2-1;
% scale yrange space of model into [-1,1]
pts(2,:) = (pts(2,:)-yrange(1))/diff(yrange)*2-1;
% scale depth space of model into [-1,1]
minz = min(pts(3,:)); maxz = max(pts(3,:));
pts(3,:) = (pts(3,:)-minz)/(maxz-minz)*2-1;


rRes = yrange(2);
cRes = xrange(2);
[rendered,z] = render_rgb_USF( pts,tri,rgb,rRes,cRes);
z = (double(z)+1)/2*(maxz-minz)+minz;

if talk
    if size(im_c,3)==1
        rendered2 = rgb2gray(rendered);
    else
        rendered2 = rendered;
    end
    figure;imshow(rendered)
    im_merge = im_c;
    im_merge(rendered2>0) = im_merge(rendered2>0)/5 + rendered2(rendered2>0)/5*4;
    figure;imshow(im_merge)
end

end

