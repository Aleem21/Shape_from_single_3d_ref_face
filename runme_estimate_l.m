clear variables
% close all
import plot3D_helper.label_axis
fusion_path = '.\data\fusion.jpg';

%% set initial variables
folder_path = '.\data\l_test\';
talk = 0;
impaths = {'band.jpg', 'ravi.jpg','david.jpg' ,'google1.jpg','baby.jpg','scarlet.jpg'};
n = numel(impaths);
f=figure;hold on
count = 1;


for i=1:n
    impath = [folder_path impaths{i}];
%% make image
[im,im_c] = read_image(impath);
%% Run face tracker
landmarks = stasm_tracker(im,talk);

if isempty(landmarks)
    continue;
end

%% Compute pose
[Rpose, Scale] = compute_pose_YaleB(landmarks, talk, im);

%% generate ref depth map
ply_path = '.\data\ref_model.ply';
% ply_path ='.\data\sphere.ply';
[dmap_ref, n_ref, N_ref] = generate_ref_depthmap(ply_path,size(im,1), size(im,2),Rpose,im);
% [dmap_ref, n_ref] = generate_ref_depthmap(ply_path,Scale, talk, 1000, 1000, Rpose,im);


%% generate ref albedo
face_database_path = '../CroppedYale';
alb_ref = get_ref_albedo(face_database_path,landmarks, size(im),fusion_path, talk)*1.5;

% alb_ref = ((get_ref_albedo(face_database_path,landmarks, size(im),fusion_path, talk)+0.5).^0.3-0.5);
% alb_ref = alb_ref*0+0.8;
alb_ref = dmap_ref*0+1;
% alb_ref = alb_ref*0+rand(size(alb_ref));

%% estimate lighting+
s = 150;
Rpose = makehgtform('scale',1/s);
l = estimate_lighting(n_ref, alb_ref, im,4,talk);
c4 = render_model('./data/sphere.ply',l,talk,Rpose,1000,1000);

l2 = estimate_lighting(n_ref, alb_ref, im,9);
c9 = render_model('./data/sphere.ply',l2,talk,Rpose,1000,1000);
% 
% l = zeros(9,1);
% l(2) = 0.3;
% c9 = render_model('./data/sphere.ply',l,talk,1000,1000);
% disp(compute_sh_cost(c9,n_ref,alb_ref,l,4))
% l2 = estimate_lighting(n_ref, alb_ref, c9,4);
% disp(compute_sh_cost(c9,n_ref,alb_ref,l2,4))
% c9new = render_model('./data/sphere.ply',l2,talk,1000,1000);
% figure;subplot(1,2,1);imagesc(c9);subplot(1,2,2);imagesc(c9new);
% disp(l')
% disp(l2')
figure(f)
subplot(2,min(n,3),count)
imshow(im_c)
subplot(2,min(n,3),count+min(n,3))
imshow(c4)
% subplot(3,n,count+2*n)
% imshow(c9)
drawnow
% depth = estimate_depth(N_ref,alb_ref,im,dmap_ref,l,50);

if mod(count,3)==0 && i<n
    f=figure;
    count = 0;
end
count = count+1;
end