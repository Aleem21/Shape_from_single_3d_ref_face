clear variables
close all
import plot3D_helper.label_axis
% depth_file = '.\data\raw_depth.png';
% depth_map = read_depthmap(depth_file,talk);


%% generate ref depth map
talk = 0;
ply_path = '..\data\ref_model.ply';
dmap_ref = generate_ref_depthmap(ply_path,talk);

%% generate ref albedo
face_database_path = './CroppedYale';
alb_ref = generate_ref_albedo(face_database_path);
