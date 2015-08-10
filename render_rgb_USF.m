function [ rendered,z ] = render_rgb_USF( pts,tri,rgb,rRes,cRes )
%RENDER_RGB_USF Summary of this function goes here
%   Detailed explanation goes here
maxz = max(pts(3,:));
minz = min(pts(3,:));
pts(3,:) = (pts(3,:)-minz)/(maxz-minz);
[r,g,b,z] = mex_render_rgb(pts,tri,rgb,cRes,rRes);
z(r==0)=nan;
rendered(:,:,1) = r';
rendered(:,:,2) = g';
rendered(:,:,3) = b';
rendered = rendered(end:-1:1,end:-1:1,:);
z = z';
z = z(end:-1:1,end:-1:1,:)*2-1;
z = (z*(maxz-minz))+minz;
% z = z*min([cRes rRes])/2;
end

