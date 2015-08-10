function make_phong_vid( s,fname )
%MAKE_PHONG_VID Summary of this function goes here
%   Detailed explanation goes here
phong.shading(s)
colormap([0.7969    0.5156    0.2617])
make_vid(s,fname);
end

