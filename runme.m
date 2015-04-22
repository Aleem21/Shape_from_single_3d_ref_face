clear variables
close all
import plot3D_helper.label_axis
% depth_file = '.\data\raw_depth.png';
% depth_map = read_depthmap(depth_file,talk);
fusion_path = '.\data\fusion.jpg';
impath = '.\data\yaleB.jpg';
%% make image
holder = rgb2gray(imread(fusion_path));
im = imread(impath);
holder(112:111+size(im,1),67:66+size(im,2)) = im;
imwrite(holder,'.\data\temp.jpg');
%% Run face tracker
[~, output] = system(['stasm4\minimal  .\data\temp.jpg .\stasm4']);
delete('.\data\temp.jpg')
landmarks = reshape(str2num(output),2,77)-repmat([67;111],1,77);

%% generate ref depth map
talk = 2;
ply_path = '..\data\ref_model.ply';
dmap_ref = generate_ref_depthmap(ply_path,talk,1000,700);

%% generate ref albedo
face_database_path = '../CroppedYale';
alb_ref = generate_ref_albedo(face_database_path);

figure;imshow(imread(impath));hold on
pts = [33 43 53 63];
plot(landmarks(1,pts),landmarks(2,pts),'o')
figure;imshow(alb_ref);