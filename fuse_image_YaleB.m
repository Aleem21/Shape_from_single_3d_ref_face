function [  fused ] = fuse_image_YaleB( fusion_path, im )

fused = im2double(rgb2gray(imread(fusion_path)));
fused(112:111+size(im,1),67:66+size(im,2)) = im;

end

