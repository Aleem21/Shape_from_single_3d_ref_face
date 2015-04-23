function [ albedo ] = get_ref_albedo( face_database_path )
%GET_DEPTHMAP Checks if average albedo already exits. If it does, returns
%in, else computes one, saves in and returns it

% adds trailing slash at end of path if one doesnt exist already
if face_database_path(end)~='\' || face_database_path(end)~='/'
    face_database_path(end+1) = '\';
end

try
    albedo = imread([face_database_path 'avg_alb.tiff']);
catch Err
    if strcmp(Err.identifier,'MATLAB:imagesci:imread:fileDoesNotExist')==1
        albedo = generate_ref_albedo(face_database_path);
        imwrite(albedo, [face_database_path 'avg_alb.tiff']);
    end
end


end

