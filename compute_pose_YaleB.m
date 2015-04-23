function [ sRT ] = compute_pose_YaleB( pp )
%COMPUTE_POSE_YALEB computes pose for specific yaleB dataset with our
%specific model. p are the landmark points in form of 2x4 matrix. Points
%being in order of left eye, right eye, nose, middle of top lip. The
%origin is on top left.

%% data conditioning
% flipping of y axis to set +y on top and hence origin at bottom left
% corner of image
sz = 192;
pp(2,:) = sz -pp(2,:);

%% 3D points for Our model
Pp = [-0.3769    0.2345  -0.08816    -0.07125;
       0.1458    0.2297  -0.05439    -0.3504;
       0.5042    0.5717   0.8754      0.866];

%% Pose computation
sRT = compute_pose(pp,Pp);

end

