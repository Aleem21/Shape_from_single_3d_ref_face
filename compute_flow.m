function [ o1_out,o2_out,o3_out,uv ] = compute_flow( im2,im1,landmarks,o1,o2,o3,n_levels,n_iters,morph,talk )
%COMPUTE_FLOW Summary of this function goes here
%   Detailed explanation goes here
if nargin<6
    talk = 0;
end


im1(isnan(im1)) = 0;
im2(isnan(im2)) = 0;

if morph
    eye_l = 31:38;
    eye_r = 41:48;

    eyes_up_i = [29 30]-16;
    eyes_center_i = [39 40]-16;
    valid = 17:77;
    valid([eyes_center_i eyes_up_i]) = [];
    
    eye_l(end)=[]; eye_r(end)=[];
    landmarks_ref = stasm_tracker(im1,4);
    p = [landmarks_ref(:,valid)...
        mean(landmarks_ref(:,eye_l),2) mean(landmarks_ref(:,eye_r),2)];
    q = [landmarks(:,valid)...
        mean(landmarks(:,eye_l),2) mean(landmarks(:,eye_r),2)];
    step = 1;
    [X,Y] = meshgrid(1:step:size(im1,2),1:step:size(im1,1));
    gv = [X(:)';Y(:)'];
    cd 'MLS'
    mlsd = MLSD2DpointsPrecompute(p,gv);
    imgo = MLSD2DWarp(im1,mlsd,q,X,Y);
    %     figure; imshow(im1);
    %     figure; imshow(im2);
    %     figure; imshow(imgo);
    o1 = MLSD2DWarp(o1,mlsd,q,X,Y);
    o2 = MLSD2DWarp(o2,mlsd,q,X,Y);
    o3 = MLSD2DWarp(o3,mlsd,q,X,Y);
    
    im1 = imgo;
    cd ..
end

if n_levels>0 && n_iters>0
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
            for j=1:max(n_iters-(n_iters/(numel(im1_pyr)-1)*(i-1)),1)
                im1_cur  = interp2(im1_pyr{i}, x-uv(:,:,1), y-uv(:,:,2));
                im1_cur(isnan(im1_cur)) = im1_pyr{i}(isnan(im1_cur));
                
                uv_cur = flow.mex_LDOF(repmat(im1_cur,1,1,3)*255,repmat(im2_pyr{i},1,1,3)*255);
                uv = uv + uv_cur;
            end
        end
        %     uv = flow.mex_LDOF(repmat(im1,1,1,3)*255,repmat(im2,1,1,3)*255);
    else
        uv = flow.mex_LDOF(im1*255,im2*255);
    end
else
    uv = zeros(size(im1,1),size(im1,2),2);
end
% gauss=fspecial('gaussian',5,5);

gauss = 1;




% generate synthetic test data, for experimenting
[x, y] = meshgrid(1:size(im1,2), 1:size(im1,1));
vx = conv2(uv(:,:,1),gauss,'same');   % an arbitrary flow field, in this case
vy = conv2(uv(:,:,2),gauss,'same');   % representing shear
uv(:,:,1) = vx; uv(:,:,2)=vy;
% compute the warped image - the subtractions are because we're specifying
% where in the original image each pixel in the new image comes from
im3 = interp2(double(im1), x-vx, y-vy);
im3(isnan(im3))=0;
im3 = real(im3);
% display the result
if talk
    %     figure;imshow(im1);
    %     figure;imshow(im2);
    err = abs(im2-im3); err(isnan(err)) = 0;
    err_val = sum(err(:));
    figure;imshow(im3);
    title(sprintf('levels=%d, iters=%d, morphing=%d, error=%d',n_levels,n_iters,morph,err_val))
    figure;imagesc(err);
    title(sprintf('levels=%d, iters=%d, morphing=%d, error=%d',n_levels,n_iters,morph,err_val))
    figure;imshow(imfuse(im2,im3,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]))
    title(sprintf('levels=%d, iters=%d, morphing=%d, error=%d',n_levels,n_iters,morph,err_val))
    
    %     cd 'black'
    %     addpath(genpath('.\utils'))
    %     figure; subplot(1,2,1);imshow(uint8(flowToColor(uv))); title('Middlebury color coding');
    %     subplot(1,2,2); plotflow(uv);   title('Vector plot');
    %     cd '..'
end


o1_out = interp2(o1, x-vx, y-vy);
% o1_out(isnan(o1_out)) = o1(isnan(o1_out));
o2_out = interp2(o2, x-vx, y-vy);
% o2_out(isnan(o2_out)) = o2(isnan(o2_out));
o3_out = interp2(o3, x-vx, y-vy);
% o3_out(isnan(o3_out)) = o3(isnan(o3_out));
end

