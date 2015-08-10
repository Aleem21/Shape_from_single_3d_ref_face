function make_phong_vid( s,fname )
%MAKE_PHONG_VID Summary of this function goes here
%   Detailed explanation goes here
phong.shading(s)
colormap([1 0.7 0.1;1 0.7 0.1])
make_vid(s,fname);
end

