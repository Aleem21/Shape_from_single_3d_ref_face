function [ depth ] = compute_depthmap( pts, tri, xrange, yrange, cRes, rRes,center )
%COMPUTE_DEPTHMAP using the pts and triangulation, generates the depthmap
%from canonical point of view in orthograhpic projection
%
%Theoratically, we are shooting rays from behind the model (from z=0 plane)
%and we check for the intersection with largest z value. This makes sure we
%get the point closest from the opposite end i.e. closest to the observer
%placed on positive z axis looking towards origin.

if nargin < 7
    center = 1;
end
if cRes<200
    scale = 200/cRes;
    cRes2 = 200;
    rRes2 = ceil(rRes*scale);
else
    cRes2 = cRes;
    rRes2 = rRes;
end
    
%% Data conditioning
if size(tri,1)>4
    tri = tri';
end
if size(pts,1)>4
    pts = pts';
end
%% conditioning
ptsn = pts;
if center
    ptsn(1,:) = (ptsn(1,:)/cRes-0.5)*2;
    ptsn(2,:) = (ptsn(2,:)/rRes-0.5)*2;
end
ptsn(3,:) = -ptsn(3,:);

tri = tri-1;

%% apply the mex function
depth  = (mex_compute_depth(ptsn, tri, xrange,yrange,cRes2, rRes2)-0.5)*2;
depth = -depth(:,end:-1:1)';
depth(depth==-1 | depth==1) = nan;
if cRes2~=cRes
    depth = imresize(depth,[rRes,cRes]);
% 
% depth  = (mex_compute_depth_new(ptsn, tri, xrange,yrange,cRes, rRes)*2)-1;
% depth = -depth(:,end:-1:1)';
% depth(depth==1) = nan;
end

