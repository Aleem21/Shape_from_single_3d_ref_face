function depth_map = read_depthmap (imname,talk)
% READ_DEPTHMAP reads the depthmap from png file and applies some processing like
% cropping,scalling.

if nargin <2
    talk = 0;
end

[~,~,a] = imread(imname,'png');
a = double(a);
a(a==0)=255;

depth_map = a(300:1450,150:1050)*1.7;

if talk
    figure;
    surf(depth_map,'edgealpha',0)
    axis equal
end