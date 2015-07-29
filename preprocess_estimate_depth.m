function [ face,face_inds, inface_inds,in_face,r_face,c_face,...
    r_inface,c_inface, b_out_full,b_in_full ] =...
    preprocess_estimate_depth( z_ref )
%PREPROCESS_ESTIMATE_DEPTH Summary of this function goes here
%   Detailed explanation goes here

face = ~isnan(z_ref);
face = remove_bad_boundary(face)>0;


% face inds
[r_face,c_face] = find(face);
face_inds = sub2ind(size(face),r_face,c_face);


% boundary
% [b_out,b_in] = find_boundary(face,false);
[b_out_full,b_in_full] = find_boundary(face,true);

% inface
in_face = (face-b_out_full)>0;
[r_inface,c_inface] = find(in_face);
inface_inds = sub2ind(size(face),r_inface,c_inface);


end

