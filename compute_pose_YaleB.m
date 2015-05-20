function [ sRT,S ] = compute_pose_YaleB( pp,talk,im )
%COMPUTE_POSE_YALEB computes pose for specific yaleB dataset with our
%specific model. p are the landmark points in form of 2x4 matrix. Points
%being in order of left eye, right eye, nose, middle of top lip. The
%origin is on top left.

if nargin<2
    talk = 0;
end

pts = [33 43 53 63 75];
pp = pp(:,pts);

%% data conditioning
% flipping of y axis to set +y on top and hence origin at bottom left
% corner of image
sz = size(im,1);
pp(2,:) = sz -pp(2,:);

%% 3D points for Our model
Pp = [-0.3769    0.2345  -0.08995    -0.06465   -0.0228;
    0.1458    0.2297  -0.03936    -0.3102    -0.4931;
    0.5042    0.5717   0.8696      0.8524     0.8576];

%% Pose computation
[sRT,S] = compute_pose(pp,Pp);


if talk
    if nargin<3
        disp('Not displaying pose matching. No input image provided');
    else
        figure;imshow(im);hold on
        plot(pp(1,:),sz - pp(2,:),'o')
        truesize
    end
    
end

end

