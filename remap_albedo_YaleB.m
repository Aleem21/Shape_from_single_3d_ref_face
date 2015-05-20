function [ albedo_trans ] = remap_albedo_YaleB( fusion_path, landmarks,albedo,sz_im )
%REMAP_ALBEDO_YALEB Summary of this function goes here
%   Detailed explanation goes here
import cv.homography_solve

fused = fuse_image_YaleB(fusion_path,albedo);
landmarks_ref = stasm_tracker_YaleB(fused,0);
if isempty(landmarks_ref)
    warning('No landmarks found in ref albedo. Using the unaltered version now.')
    return
end
pts = [33 43 53 63 75];
% pts = 1:size(landmarks_ref,2);
H = homography_solve(landmarks_ref(:,pts),landmarks(:,pts));
H = H';
% H(:,1:2) = H(:,[2 1]);
% H(1:2,:) = H([2 1],:);
Ht = maketform('Projective',H);
albedo_trans = imtransform(albedo,Ht,'XData',[1 sz_im(2)],'YData',[1 sz_im(1)]);
albedo_trans = imresize(albedo_trans,sz_im);
end

