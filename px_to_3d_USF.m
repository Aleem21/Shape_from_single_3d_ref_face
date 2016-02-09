function [ pts3D ] = px_to_3d_USF( pts2D,depthmap,cRes,rRes )
%PX_TO_3D_USF Summary of this function goes here
%   Detailed explanation goes here
% pts has to be two(rows) by number of points(cols)
if size(pts2D,1)>2
    pts2D = pts2D';
end

pts3D = -(pts2D./repmat([cRes rRes]',1,size(pts2D,2)) - 0.5) * 2;
pts3D(3,:) = depthmap(sub2ind(size(depthmap),pts2D(2,:),pts2D(1,:)))*2/(min([cRes rRes])-1);
end

