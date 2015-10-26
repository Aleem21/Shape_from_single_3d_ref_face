function [ xy ] = rgb2xy( rgb )
%CHROMATICITY_RGB2XY Converts RGB to CIE xy space going through
%RGB->XYZ(CIE)->xy(CIE)
% rgb needs to be a 3 vector

XYZ = reshape(rgb2xyz(reshape(rgb,1,1,3)),3,1,1);

xy = XYZ(1:2)/sum(XYZ);

end

