function [ sRT,S,sT ] = compute_pose_USF( pp,talk,im,restrictive,cRes,rRes)

if nargin<2
    talk = 0;
end
if nargin < 4
    restrictive = 0;
end
[pts,tri,rgb] = read_USF_eko('D:\Drives\Google Drive\Research UCSD\Ravi\Sony SFS\datasets\USF 3D Face Data\USF Raw 3D Face Data Set\data_files\test',512,512);

[rendered,z] = render_rgb_USF(pts,tri,rgb,rRes,cRes);
z = double(z)*min(rRes,cRes)/2;
landmarks_ref = stasm_tracker(rendered,talk);
valid = 17:77;
landmarks_3D_ref = px_to_3d_USF(landmarks_ref(:,valid),z,cRes,rRes);
pp = pp(:,valid);

%% data conditioning
% flipping of y axis to set +y on top and hence origin at bottom left
% corner of image
sz1 = size(im,2);
pp(1,:) = sz1 -pp(1,:);
sz2 = size(im,1);
pp(2,:) = sz2 -pp(2,:);

%% 3D points for Our model
% Pp = [-0.3769    0.2345  -0.08995    -0.06465   -0.0228;
%     0.1458    0.2297  -0.03936    -0.3102    -0.4931;
%     0.5042    0.5717   0.8696      0.8524     0.8576];
Pp = landmarks_3D_ref;
Pp(1,:) = Pp(1,:);
Pp(2,:) = Pp(2,:);
%% Pose computation
[sRT,S,sT] = compute_pose(pp,Pp);
if restrictive
    sRT = sT;
end
if talk
    if nargin<3
        disp('Not displaying pose matching. No input image provided');
    else
        figure;imshow(im);hold on
        plot(sz1-pp(1,:),sz2 - pp(2,:),'ro','markersize',10)
        truesize
        
        pts2 = sRT*[pts; ones(1,size(pts,2))];
        pts2 = pts2(1:3,:);
        
        pts2(3,:) = pts(3,:)*min(S);
        figure;imshow(im); hold on
%         plot(sz1-pts2(1,:),sz2-pts2(2,:),'.')
        landmarks2 = sRT*[landmarks_3D_ref ;ones(1,size(landmarks_3D_ref,2))];
        plot(sz1-landmarks2(1,:),sz2-landmarks2(2,:),'or')
        plot(sz1-pp(1,:),sz2-pp(2,:),'.g')
        title('green: Target landmarks, red: reprojected 3D model landmarks')
        
    end
    
end





end

