function [ pts_rotated,tri,clr ] = model2GL( model,midx,midy )
%MODEL2GL Takes arguments from model fitting function and converts them
%into a form suitable for openGL depth finding function
% midx is mean of xrange usually

x2 = model.R*model.FV.vertices(:,:)';
x2(1,:) = x2(1,:)+model.t(1);
x2(2,:)=x2(2,:)+model.t(2);
x2 = x2*model.s;

pts_rotated = double(x2);
pts_rotated(1,:) = -(pts_rotated(1,:) - midx) + midx;
%pts_rotated(2,:) = -(pts_rotated(2,:) - midy) + midy;

pts_rotated(3,:) = -(pts_rotated(3,:) - max(pts_rotated(3,:)));

tri = model.FV.faces'-1;
if nargout>2
    load('mask')
    clr = double(repmat(mask',3,1)==0);
end
end

