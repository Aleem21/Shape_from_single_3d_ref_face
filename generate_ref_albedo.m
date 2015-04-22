function [ avg_albedo ] = generate_ref_albedo( face_database_path )
%GENERATE_REF_ALBEDO Summary of this function goes here
%   Detailed explanation goes here

impaths = cell(1,0);
for i=1:39
    if i<10
        no = ['0' num2str(i)];
    else
        no = num2str(i);
    end
    impaths{i} = [face_database_path '/yaleB' no '/yaleB' no '_P00A+000E+00.pgm'];    
end

landmarks = cell(1,numel(impaths));


if numel(impaths)>0
    landmarks_dst = compute_features(impaths{15});
    im1 = im2double(imread(impaths{15}));
    avg_albedo = im1;
end

sum_map = ones(size(im1));

for i=2:numel(impaths)
    landmarks_src = compute_features(impaths{i});
    im_i = im2double(imread(impaths{i}));
    im_fixed = rescale_face(landmarks_dst,landmarks_src,im_i);
    avg_albedo = avg_albedo + im_fixed;
    sum_map = sum_map + (im_fixed>0);
end
avg_albedo = avg_albedo./sum_map;

end

