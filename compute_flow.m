function [ o1_out,o2_out,o3_out,uv ] = compute_flow( im2,im1,landmarks,o1,o2,o3,n_levels,n_iters,morph,talk,yale_landmarks )
%COMPUTE_FLOW Summary of this function goes here
%   Detailed explanation goes here
if nargin<10
    talk = 0;
end
if nargin<11
    yale_landmarks = [];
end


im1(isnan(im1)) = 0;
im2(isnan(im2)) = 0;

if morph
    addpath('./interparc/')
    eye_l = 31:38;
    eye_r = 41:48;

    eyes_up_i = [29 30]-16;
    eyes_center_i = [39 40]-16;
    eye_b_l = 17:22;
    eye_b_r = 23:28;
    
    valid = 17:77;
    valid([eyes_center_i eyes_up_i]) = [];
    valid(ismember(valid,eye_l) | ismember(valid,eye_r)) = [];
    valid(ismember(valid,eye_b_l) | ismember(valid,eye_b_r)) = [];
    l_u = valid(12:21);
    l_l = valid([12 22:24 18 25:29]);
    valid(ismember(valid,l_u) | ismember(valid,l_l)) = [];

    eye_l(end)=[]; eye_r(end)=[];
    if ~isempty(yale_landmarks)
        landmarks_ref = yale_landmarks;
    else
        landmarks_ref = stasm_tracker(imresize(im1,2),4);
        landmarks_ref = floor(landmarks_ref/2);
    end
    %% eye brows
    if numel(landmarks_ref)>0
        p_eye_b_l_pt = landmarks_ref(:,eye_b_l([6 1 2 3]))/2 + landmarks_ref(:,eye_b_l([6 5 4 3]))/2;
        p_eye_b_r_pt = landmarks_ref(:,eye_b_r([1 2 3 4]))/2 + landmarks_ref(:,eye_b_r([1 6 5 4]))/2;
        p_eye_b_l_pt = interparc([0:0.01:1],p_eye_b_l_pt(1,:),p_eye_b_l_pt(2,:))';
        p_eye_b_r_pt = interparc([0:0.01:1],p_eye_b_r_pt(1,:),p_eye_b_r_pt(2,:))';

        q_eye_b_l_pt = landmarks(:,eye_b_l([6 1 2 3]))/2 + landmarks(:,eye_b_l([6 5 4 3]))/2;
        q_eye_b_r_pt = landmarks(:,eye_b_r([1 2 3 4]))/2 + landmarks(:,eye_b_r([1 6 5 4]))/2;
        q_eye_b_l_pt = interparc([0:0.01:1],q_eye_b_l_pt(1,:),q_eye_b_l_pt(2,:))';
        q_eye_b_r_pt = interparc([0:0.01:1],q_eye_b_r_pt(1,:),q_eye_b_r_pt(2,:))';

        %% lips
        p_l_u = landmarks_ref(:,l_u);
%         sp_l_u = [interparc([0:0.02:1],p_l_u(1,[1:7]),p_l_u(2,[1:7]))'...
%             interparc([0:0.02:1],p_l_u(1,[7:end 1]),p_l_u(2,[7:end 1]))'];
        sp_l_u = [interparc([0:0.02:1],p_l_u(1,[1:7]),p_l_u(2,[1:7]))'...
            ];

        p_l_l = landmarks_ref(:,l_l);
%         sp_l_l = [interparc([0:0.02:1],p_l_l(1,[1:5]),p_l_l(2,[1:5]))'...
%             interparc([0:0.02:1],p_l_l(1,[5:end 1]),p_l_l(2,[5:end 1]))'];
        sp_l_l = [interparc([0:0.02:1],p_l_l(1,[5:10 1]),p_l_l(2,[5:10 1]))'...
            ];
        p_lip = [sp_l_u sp_l_l];

        q_l_u = landmarks(:,l_u);
%         sq_l_u = [interparc([0:0.02:1],q_l_u(1,[1:7]),q_l_u(2,[1:7]))'...
%             interparc([0:0.02:1],q_l_u(1,[7:end 1]),q_l_u(2,[7:end 1]))'];
        sq_l_u = [interparc([0:0.02:1],q_l_u(1,[1:7]),q_l_u(2,[1:7]))'...
            ];
        q_l_l = landmarks(:,l_l);
%         sq_l_l = [interparc([0:0.02:1],q_l_l(1,[1:5]),q_l_l(2,[1:5]))'...
%             interparc([0:0.02:1],q_l_l(1,[5:end 1]),q_l_l(2,[5:end 1]))'];
        sq_l_l = [interparc([0:0.02:1],q_l_l(1,[5:10 1]),q_l_l(2,[5:10 1]))'...
            ];
        q_lip = [sq_l_u sq_l_l];

        %% eyes
        p_e_l = landmarks_ref(:,eye_l);
    %     sp_e_l = [interparc([0:0.02:1],p_e_l(1,[1:5]),p_e_l(2,[1:5]))'...
    %         interparc([0:0.02:1],p_e_l(1,[5:end 1]),p_e_l(2,[5:end 1]))'];
        sp_e_l = [interparc([0:0.02:1],p_e_l(1,[1:5]),p_e_l(2,[1:5]))'];

        p_e_r = landmarks_ref(:,eye_r);
    %     sp_e_r = [interparc([0:0.02:1],p_e_r(1,[1:5]),p_e_r(2,[1:5]))'...
    %         interparc([0:0.02:1],p_e_r(1,[5:end 1]),p_e_r(2,[5:end 1]))'];
        sp_e_r = [interparc([0:0.02:1],p_e_r(1,[1:5]),p_e_r(2,[1:5]))'];
        p_e = [sp_e_l sp_e_r];

        q_e_l = landmarks(:,eye_l);
    %     sq_e_l = [interparc([0:0.02:1],q_e_l(1,[1:5]),q_e_l(2,[1:5]))'...
    %         interparc([0:0.02:1],q_e_l(1,[5:end 1]),q_e_l(2,[5:end 1]))'];
        sq_e_l = [interparc([0:0.02:1],q_e_l(1,[1:5]),q_e_l(2,[1:5]))'];

        q_e_r = landmarks(:,eye_r);
    %     sq_e_r = [interparc([0:0.02:1],q_e_r(1,[1:5]),q_e_r(2,[1:5]))'...
    %         interparc([0:0.02:1],q_e_r(1,[5:end 1]),q_e_r(2,[5:end 1]))'];
        sq_e_r = [interparc([0:0.02:1],q_e_r(1,[1:5]),q_e_r(2,[1:5]))'];

        q_e = [sq_e_l sq_e_r];
    %% nose 
        p = [mean(landmarks_ref(:,valid),2)...
            mean(landmarks_ref(:,eye_l),2) mean(landmarks_ref(:,eye_r),2)...
            p_eye_b_l_pt p_eye_b_r_pt p_lip p_e];
        q = [mean(landmarks(:,valid),2)...
            mean(landmarks(:,eye_l),2) mean(landmarks(:,eye_r),2)...
            q_eye_b_l_pt q_eye_b_r_pt q_lip q_e];
        if norm(p(:)-q(:))<500
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
    end
    
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
                
                uv_cur = flow.mex_LDOF(repmat(im1_cur,1,1,3)*255,repmat(im2_pyr{i},1,1,3)*255,0.8,30);
%                 uv_cur = flow.mex_LDOF(repmat(im1_cur,1,1,3)*255,repmat(im2_pyr{i},1,1,3)*255);
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

