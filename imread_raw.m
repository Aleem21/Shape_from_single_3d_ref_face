function [ lin_rgb,lin_gray ] = imread_raw( impath )
%READ_RAW Summary of this function goes here
%   Detailed explanation goes here
% system(['dcraw-9.26-ms-64-bit.exe -o 2 -w -W -6 -4 -T' impath])
% tiff_path = [impath(1:end-3) 'tiff'];
% lin_srgb = double(imread(tiff_path));
% return;
RGB_to_XYZ = [  0.4124564 0.3575761 0.1804375;
                0.2126729 0.7151522 0.0721750;
                0.0193339 0.1191920 0.9503041];
XYZ_to_cam = [   6847,-614, -1014;
                -4669, 12737,2139;
                -1197, 2488, 6846]/10000;
system (['dcraw-9.26-ms-64-bit -4 -T -W -w -q 3 -o 5 ' impath]);
tiff_path = [impath(1:end-3) 'tiff'];
lin_XYZ = im2double(imread(tiff_path));

% lin_gray = rgb2gray(xyz2rgb(lin_XYZ));
lin_gray = lin_XYZ(:,:,2);
XYZ_to_RGB = inv(RGB_to_XYZ);
lin_rgb = correct(lin_XYZ,XYZ_to_RGB);

return;
system (['dcraw-9.26-ms-64-bit -4 -D -T ' impath]);
tiff_path = [impath(1:end-3) 'tiff'];
raw = double(imread(tiff_path));
try
    meta_info  = imfinfo([impath(1:end-3) 'dng']);
    black = meta_info.SubIFDs{1}.BlackLevel(1);
    saturation = meta_info.SubIFDs{1}.WhiteLevel;
    wb_multipliers = (meta_info.AsShotNeutral).^-1;
    XYZ_to_cam = reshape(meta_info.ColorMatrix2,3,3);
    XYZ_to_cam = XYZ_to_cam';
catch
    [~,header] = system (['dcraw-9.26-ms-64-bit -v -W -T ' impath]);
    d_i = strfind(header,'darkness ');
    d_end =  strfind(header(d_i+9:end),',');
    black = str2double(header(d_i+8 +(1:d_end(1)-1)));
    s_i = d_i+8+d_end(1)+13;
    s_end = strfind(header(s_i:end),' ');
    saturation = str2double(header(s_i + (0:s_end(1)-3)));
    wb_i = strfind(header,'multipliers');
    endl = 10;
    wb_end = strfind(header(wb_i:end),endl);
    wb_multipliers = str2num(header(wb_i + (12:wb_end-2)));
    wb_multipliers = wb_multipliers(1:3);
    
    
end
lin_bayer = (raw-black)/(saturation-black);
lin_bayer = max(0,min(lin_bayer,1));
mask = wbmask(size(lin_bayer,1),size(lin_bayer,2),wb_multipliers,'gbrg');
balanced_bayer  = lin_bayer.*mask;

temp = uint16(balanced_bayer/max(balanced_bayer(:))*2^16);
lin_cam = double(demosaic(temp,'gbrg'))/2^16;


RGB_to_cam = XYZ_to_cam*RGB_to_XYZ;
RGB_to_cam = RGB_to_cam./ repmat(sum(RGB_to_cam,2),1,3);
cam_to_RGB = inv(RGB_to_cam);
lin_rgb = correct(lin_cam,cam_to_RGB);

cam_to_XYZ = inv(XYZ_to_cam);
lin_XYZ = correct(lin_cam,cam_to_XYZ);
lin_gray = lin_XYZ(:,:,2);


end
function colormask = wbmask(m,n,wbmults,align)
% COLORMASK = wbmask(M,N,WBMULTS,ALIGN)
%
% Makes a white-balance multiplicative mask for an image of size m-by-n
% with RGB while balance multipliers WBMULTS = [R_scale G_scale B_scale].
% ALIGN is string indicating Bayer arrangement: ’rggb’,’gbrg’,’grbg’,’bggr’
colormask = wbmults(2)*ones(m,n); %Initialize to all green values
switch align
    case 'rggb'
        colormask(1:2:end,1:2:end) = wbmults(1); %r
        colormask(2:2:end,2:2:end) = wbmults(3); %b
    case 'bggr'
        colormask(2:2:end,2:2:end) = wbmults(1); %r
        colormask(1:2:end,1:2:end) = wbmults(3); %b
    case 'grbg'
        colormask(1:2:end,2:2:end) = wbmults(1); %r
        colormask(1:2:end,2:2:end) = wbmults(3); %b
    case 'gbrg'
        colormask(2:2:end,1:2:end) = wbmults(1); %r
        colormask(1:2:end,2:2:end) = wbmults(3); %b
end
end
function corrected = correct(I,correction)
sz = size(I);
corrected = reshape((correction*(reshape(I,[],3)'))',sz);
corrected = max(0,min(corrected,1));
end
