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

if numel(impaths)>0
    im1 = im2double(imread(impaths{1}));
    avg_albedo = im1;
end

for i=2:numel(impaths)
    im_i = im2double(imread(impaths{i}));
    avg_albedo = avg_albedo*(i-1)/i + im_i/i;
end

end

