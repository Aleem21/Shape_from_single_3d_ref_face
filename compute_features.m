function [ landmarks ] = compute_features( impath,talk )
%COMPUTE_FEATURES puts im_path image in a placeholder to compute complete
%face and runs facetracker stasm4.0 to compute locations of eyes, nose and
%middle of top lip
if nargin<2
    talk = 0;
end

fusion_path = '.\data\fusion.jpg';
loc = [112 67];
%% make and write temp(fullface) image for face tracker
im = imread(impath);
if size(im,3)==3
    im = rgb2gray(im);
end
holder = rgb2gray(imread(fusion_path));
holder([0:size(im,1)-1]+loc(1),[0:size(im,2)-1]+loc(2)) = im;
imwrite(holder,'.\data\temp.jpg');

%% Run face tracker
[~, output] = system(['stasm4\minimal  .\data\temp.jpg .\stasm4']);
delete('.\data\temp.jpg')

%% extract features
try
    landmarks = reshape(str2num(output),2,77)-repmat([67;111],1,77);
    pts = [33 43 53 63];
    landmarks = landmarks(:,pts);
catch
    landmarks = [];
    return;
end
%% display results
if talk
    imshow(imread(impath));hold on
    plot(landmarks(1,:),landmarks(2,:),'o')
end

end

