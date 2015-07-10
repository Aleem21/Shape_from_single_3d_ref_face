function [ depth ] = compute_general_depthmap( pts, tri, cRes, rRes, xrange, yrange)
%COMPUTE_DEPTHMAP Computes depth map using compute_depthmap_unit working as
% a wrapper function, and giving minimum depth equal to 1 and keeping the
% aspect ratio to be equal to 1 for x,y,z axis.

aspectRatio_model = diff(xrange)/diff(yrange);
aspectRatio_img = cRes/rRes;

if (aspectRatio_model~=aspectRatio_img)
    warning(['Aspect ratio of model range and image is not same. ' ...
    'The output depth map aspect ratio might be scewed.']);
end

ratio_model_to_img = cRes/diff(xrange);

pts_scaled = pts;


% scale xrange space of model into [-1,1]
pts_scaled(1,:) = (pts_scaled(1,:)-xrange(1))/diff(xrange)*2-1;
% scale yrange space of model into [-1,1]
pts_scaled(2,:) = (pts_scaled(2,:)-yrange(1))/diff(yrange)*2-1;
% scale depth space of model into [-1,1]
minz = min(pts_scaled(3,:)); maxz = max(pts_scaled(3,:));
pts_scaled(3,:) = (pts_scaled(3,:)-minz)/(maxz-minz)*2-1;



res = max([rRes cRes]);

depth = compute_depthmap_unit(pts_scaled,tri,res,res);

%rescale depth map to match required cRes rRes
depth = imresize(depth,[rRes cRes]);

%rescale depth according to image scale
depth = depth * (maxz-minz) * ratio_model_to_img;
