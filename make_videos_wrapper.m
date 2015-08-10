function make_videos_wrapper( depth,depth_gt,im,f )
%MAKE_VIDEOS Summary of this function goes here
%   Detailed explanation goes here
if nargin<4
    f=figure;
    drawnow
    % jFrame = get(handle(gcf),'JavaFrame');
    % jFrame.setMaximized(true);
end

figure(f)
reset(f)
depth_s = surf(depth,im,'edgealpha',0,'facecolor','interp');axis equal
colormap 'gray';
phong.shading(depth_s);
make_phong_vid(depth_s,'vid_skin.mp4');

figure(f)
reset(f)
depth_s = surf(depth,im,'edgealpha',0,'facecolor','interp');axis equal
colormap 'gray';
phong.shading(depth_s);
drawnow
make_vid(depth_s,'vid_gray.mp4')


figure(f)
reset(f)
depth_s = surf(depth_gt,im,'edgealpha',0,'facecolor','interp');axis equal
colormap 'gray';
phong.shading(depth_s);
drawnow
make_phong_vid(depth_s,'vid_skin_gt.mp4')

figure(f)
reset(f)
depth_s = surf(depth_gt,im,'edgealpha',0,'facecolor','interp');axis equal
colormap 'gray';
phong.shading(depth_s);
drawnow
make_vid(depth_s,'vid_grey_gt.mp4')

end

