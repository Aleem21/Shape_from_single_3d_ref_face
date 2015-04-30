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

%% Compute pose
pts = [33 43 53 63 75];
landmarks = landmarks(:,pts);
talk = 0;
[Rpose, Scale] = compute_pose_YaleB(landmarks, talk, im);

%% generate ref depth map
talk = 2;
ply_path = '..\data\ref_model.ply';
[dmap_ref, n_ref] = generate_ref_depthmap(ply_path,Scale, talk, size(im,1), size(im,2), Rpose,im);

%% generate ref albedo
face_database_path = '../CroppedYale';
talk = 0;
alb_ref = get_ref_albedo(face_database_path, talk);


