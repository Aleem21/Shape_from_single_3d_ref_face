function [ o1_out,o2_out,o3_out,uv ] = compute_flow( im2,im1,o1,o2,o3,talk )
%COMPUTE_FLOW Summary of this function goes here
%   Detailed explanation goes here
if nargin<6
    talk = 0;
end
n_levels = 10;
im1(isnan(im1)) = 0;
im2(isnan(im2)) = 0;
cd 'black'
addpath(genpath('.\utils'))
% method = 'classic+nl-fast';
% uv = estimate_flow_interface(im1*255, im2*255, method, []);

        f1=figure;
        f2=figure;
        f3=figure;
if size(im1,3)==1
    im1_pyr{1} = im1;
    im2_pyr{1} = im2;
    for i=2:n_levels
        if min([size(im1_pyr{i-1},1) size(im1_pyr{i-1},2)])<30
            break
        end
        im1_pyr{i} = impyramid(im1_pyr{i-1},'reduce');
        im2_pyr{i} = impyramid(im2_pyr{i-1},'reduce');
    end
    uv = zeros(size(im1_pyr{end},1),size(im1_pyr{end},2),2);
    uv_cur = zeros(size(im1_pyr{end},1),size(im1_pyr{end},2),2);
    im1_cur = im1_pyr{end};
    for i=numel(im1_pyr):-1:1
        uv = imresize(uv,[size(im1_pyr{i},1) size(im1_pyr{i},2)]);
        [x, y] = meshgrid(1:size(im1_pyr{i},2), 1:size(im1_pyr{i},1));

        im1_cur  = interp2(im1_pyr{i}, x-uv(:,:,1), y-uv(:,:,2));
        im1_cur(isnan(im1_cur)) = im1_pyr{i}(isnan(im1_cur));
        figure(f1);imshow(im1_pyr{i})
        figure(f2);imshow(im2_pyr{i})
        figure(f3);imshow(im1_cur)

        uv_cur = flow.mex_LDOF(repmat(im1_cur,1,1,3)*255,repmat(im2_pyr{i},1,1,3)*255);
         uv = uv + uv_cur;
    end
%     uv = flow.mex_LDOF(repmat(im1,1,1,3)*255,repmat(im2,1,1,3)*255);
else
    uv = flow.mex_LDOF(im1*255,im2*255);
end
gauss=fspecial('gaussian',5,5);

gauss = 1;




% generate synthetic test data, for experimenting
[x, y] = meshgrid(1:size(im1,2), 1:size(im1,1));
vx = conv2(uv(:,:,1),gauss,'same');   % an arbitrary flow field, in this case
vy = conv2(uv(:,:,2),gauss,'same');   % representing shear
uv(:,:,1) = vx; uv(:,:,2)=vy;
% compute the warped image - the subtractions are because we're specifying
% where in the original image each pixel in the new image comes from
im3 = interp2(double(im1), x-vx, y-vy);

% display the result
if talk
    figure;imshow(im1)
    figure;imshow(im2)
    figure;
    imshow(im3);
    figure; subplot(1,2,1);imshow(uint8(flowToColor(uv))); title('Middlebury color coding');
    subplot(1,2,2); plotflow(uv);   title('Vector plot');
    
end

cd '..'

o1_out = interp2(o1, x-vx, y-vy);
% o1_out(isnan(o1_out)) = o1(isnan(o1_out));
o2_out = interp2(o2, x-vx, y-vy);
% o2_out(isnan(o2_out)) = o2(isnan(o2_out));
o3_out = interp2(o3, x-vx, y-vy);
% o3_out(isnan(o3_out)) = o3(isnan(o3_out));
end

