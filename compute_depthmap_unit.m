function [ depth ] = compute_depthmap_unit( pts, tri, rRes, cRes)
%COMPUTE_DEPTHMAP using the pts and triangulation, generates the depthmap
%from canonical point of view in orthograhpic projection
%
%Theoratically, we are shooting rays from behind the model (from z=0 plane)
%and we check for the intersection with largest z value. This makes sure we
%get the point closest from the opposite end i.e. closest to the observer
%placed on positive z axis looking towards origin.
if nargin < 4
    rRes = 400;
end
if nargin < 3
    cRes = 400;
end

%% initializations and conditioning

% legacy variables
xrange = [-1 1];
yrange = [-1 1];

% if resolution is smaller than this specific value, openGL behaves weirdly
if cRes<200
    scale = 200/cRes;
    cRes2 = 200;
    rRes2 = ceil(rRes*scale);
else
    cRes2 = cRes;
    rRes2 = rRes;
end

% Data conditioning
if size(tri,1)>4
    tri = tri';
end
if size(pts,1)>4
    pts = pts';
end

pts(3,:) = -pts(3,:);
tri = tri-1;

%% apply the mex function to compute depth through a call to OPEN GL

depth  = (mex_compute_depth(pts, tri, xrange,yrange,cRes2, rRes2)-0.5)*2;
depth = -depth(:,end:-1:1)';
depth(depth==-1 | depth==1) = nan;
if cRes2~=cRes
    depth = imresize(depth,[rRes,cRes]);
end

% convert depth into the inverse map, i.e. the lesser the value, the
% closer it is from you (camera, placed at origin)
depth = 1-depth;
% depth = 0 meanns corresponds to z=1 in the model, depth = 2 corresponds
% to z=-1 in model
