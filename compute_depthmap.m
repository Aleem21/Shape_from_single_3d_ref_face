function [ depth ] = compute_depthmap( pts, tri, xrange, yrange, cRes, rRes )
%COMPUTE_DEPTHMAP using the pts and triangulation, generates the depthmap
%from canonical point of view in orthograhpic projection
%
%Theoratically, we are shooting rays from behind the model (from z=0 plane)
%and we check for the intersection with largest z value. This makes sure we
%get the point closest from the opposite end i.e. closest to the observer
%placed on positive z axis looking towards origin.

%% Data conditioning
if size(tri,1)>4
    tri = tri';
end
if size(pts,1)>4
    pts = pts';
end

%% conditioning
ptsn = pts;
ptsn(1,:) = (ptsn(1,:)/cRes-0.5)*2;
ptsn(2,:) = (ptsn(2,:)/rRes-0.5)*2;
ptsn(3,:) = -ptsn(3,:);

tri = tri-1;

%% apply the mex function
depth = (mex_compute_depth(ptsn, tri, [-1 1],[-1 1],cRes, rRes)-0.5)*2;
depth = 1-depth(:,end:-1:1)';
depth(depth==0) = nan;
end

