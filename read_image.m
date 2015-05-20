function [ im,im_c ] = read_image( impath )
%READ_IMAGE Summary of this function goes here
%   Detailed explanation goes here


im_c = im2double(imread(impath));
if size(im_c,3)==3
    im = rgb2gray(im_c);
else
    im = im_c;
end




end

