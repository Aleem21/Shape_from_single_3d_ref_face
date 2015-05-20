function [ albedo ] = get_ref_albedo_YaleB( face_database_path, landmarks,sz_im,fusion_path,talk )
%GET_DEPTHMAP Checks if average albedo already exits. If it does, returns
%in, else computes one, saves in and returns it
if nargin<4
    talk = 0;
end
% adds trailing slash at end of path if one doesnt exist already
if face_database_path(end)~='\' || face_database_path(end)~='/'
    face_database_path(end+1) = '\';
end

try
    albedo = im2double(imread([face_database_path 'avg_alb.tiff']));
catch Err
    if strcmp(Err.identifier,'MATLAB:imagesci:imread:fileDoesNotExist')==1
        albedo = generate_ref_albedo_YaleB(face_database_path);
        imwrite(albedo, [face_database_path 'avg_alb.tiff']);
    end
end
% albedo = remap_albedo_YaleB( fusion_path, landmarks,albedo,sz_im );
albedo = albedo/max(albedo(:));
if talk
    figure;imshow(albedo);
end
end

