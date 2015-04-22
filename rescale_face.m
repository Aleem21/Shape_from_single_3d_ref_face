function [ im_dst ] = rescale_face( landmarks_dst, landmarks_src, im_src )
%RESCALE_FACE rescales, translates im_src such that landmarks in src match
%that of destination

eye_w_src = landmarks_src(1,2)-landmarks_src(1,1);
eye_w_dst = landmarks_dst(1,2)-landmarks_dst(1,1);

% horizontal scaling needed
xscale = eye_w_dst/eye_w_src;
n_cols = floor(size(im_src,2)*xscale);

eye2lip_src = landmarks_src(2,4) - landmarks_src(2,1);
eye2lip_dst = landmarks_dst(2,4) - landmarks_dst(2,1);

%vertical scale needed
yscale = eye2lip_dst/eye2lip_src;
n_rows = floor(size(im_src,1)*yscale);

%find left eye offset
offset = landmarks_dst(:,1)' - landmarks_src(:,1)';

im_src_scaled = imresize(im_src,[n_rows,n_cols]);
t = maketform('affine',[eye(2);offset]);

im_dst = imtransform(im_src_scaled,t,'XData',[1 size(im_src,2)],'YData',[1 size(im_src,1)]);

%% fix rotation
theta_src = atan2(landmarks_src(2,2)-landmarks_src(2,1),landmarks_src(1,2)-landmarks_src(1,1));
theta_dst = atan2(landmarks_dst(2,2)-landmarks_dst(2,1),landmarks_dst(1,2)-landmarks_dst(1,1));
theta_diff = theta_dst - theta_src;

% from left eye in dst to origin
T = eye(3);
T(1:2,3) = [landmarks_dst([2 1],1)];

% from origin to left eye in dst
Tinv = T;
Tinv(1:2,3) = -Tinv(1:2,3);

Rot = eye(3);
Rot(1:2,1:2) = [cos(theta_diff) -sin(theta_diff);
                sin(theta_diff)  cos(theta_diff)];
t = maketform('projective',(Tinv*Rot*T)');
im_dst = imtransform(im_dst,t,'XData',[1 size(im_src,2)],'YData',[1 size(im_src,1)]);
end

