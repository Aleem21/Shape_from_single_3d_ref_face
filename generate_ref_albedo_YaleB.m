function [ avg_albedo ] = generate_ref_albedo_YaleB( face_database_path )
%GENERATE_REF_ALBEDO Summary of this function goes here
%   Detailed explanation goes here
import batch.getAllFiles


% if face_database_path(end)~='\' || face_database_path(end)~='/'
%     face_database_path(end+1) = '\';
% end
% 
% impaths = cell(1,0);
% for i=1:39
%     if i<10
%         no = ['0' num2str(i)];
%     else
%         no = num2str(i);
%     end
%     impaths{i} = [face_database_path 'yaleB' no '/yaleB' no '_P00A+000E+00.pgm'];
% end




if face_database_path(end)~='\' && face_database_path(end)~='/'
    face_database_path(end+1) = '\';
end


impaths= getAllFiles(face_database_path,1);

strt = 1;
finding = 1;
if numel(impaths)>0
    while(finding)
        try
            landmarks_dst = compute_features(impaths{strt});
            im1 = im2double(imread(impaths{15}));
            avg_albedo = im1;
            finding = 0;
        catch
            strt = strt+1;
        end
    end
end

sum_map = ones(size(im1));
fprintf('......');
for i=(strt+1):numel(impaths)
    fprintf('\b\b\b\b\b\b%.2d%%...',ceil(i/numel(impaths)*100))
    try
        landmarks_src = compute_features(impaths{i});
        im_i = im2double(imread(impaths{i}));
        im_fixed = rescale_face(landmarks_dst,landmarks_src,im_i);
        avg_albedo = avg_albedo + im_fixed;
%         sum_map = sum_map + (im_fixed>0);
        sum_map = sum_map + (im_fixed>0);
    catch
    end
end
avg_albedo = avg_albedo./sum_map;

end


