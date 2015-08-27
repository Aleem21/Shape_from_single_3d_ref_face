function do_relighting( im1,l1,n,frameRate,speed,span )
%DO_RELIGHTING Summary of this function goes here
%   Detailed explanation goes here

    a0 = pi;
    a1 = 2*pi/sqrt(3);
    a2= 2*pi/sqrt(8);
    
    c0 =    sqrt(1  /(4*pi));
    c1 =    sqrt(3  /(4*pi));
    c2 = 3* sqrt(5  /(12*pi));
    
    nx = n(:,:,1);
    nx = nx(:)';        %change to row vector
    ny = n(:,:,2);
    ny = ny(:)';
    nz = n(:,:,3);
    nz = nz(:)';
    Y = [   a0*c0*     ones(size(nx));
        a1*c1*     nx;
        a1*c1*     ny;
        a1*c1*     nz;
        a2*c2*     nx.*ny;
        a2*c2*     nx.*nz;
        a2*c2*     ny.*nz;
        a2*c2/2*    (nx.^2 - ny.^2);
        a2*c2/(12)^.5*   (3*nz.^2-1) ];
    
    
    
    l1 = l1./ [a0*c0; a1*c1; a1*c1; a1*c1;  a2*c2; a2*c2; a2*c2; a2*c2/2; a2*c2/(12)^.5 ];
    l1 = l1(1:4);
    Y = Y(1:4,:);
    im1r = im1(:,:,1);
    im1g = im1(:,:,2);
    im1b = im1(:,:,3);
	
    frames = cell(1,0);
    writerObj = VideoWriter('relighting.mp4','MPEG-4');
    writerObj.FrameRate = frameRate;
    open(writerObj);
    f=figure;
%     q1 = quaternion.angle2quat(0,0,0);
    q1 = quaternion.angle2quat(0,pi/4*span,0);
%     for i=linspace(0,1,speed)
%         q_cur = quaternion.slerp(q1,q2,i,eps);
%         [z,y,x]=quaternion.quat2angle(q_cur);
%         l2 = spherical_harmonics.rotate(l1,...
%             'xrotate',x,'yrotate',y,'zrotate',z);
%         im2Vr = max(im1r(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
%         im2Vg = max(im1g(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
%         im2Vb = max(im1b(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
%         im2(:,:,1) = reshape(im2Vr,[size(im1,1) size(im1,2)] );
%         im2(:,:,2) = reshape(im2Vg,[size(im1,1) size(im1,2)] );
%         im2(:,:,3) = reshape(im2Vb,[size(im1,1) size(im1,2)] );
%         imshow(im2);
%         figure(f)
%         drawnow
%         writeVideo(writerObj,getframe);
%     end
%     q1 = q2;
    q2 = quaternion.angle2quat(0,pi/4*span,pi/4*span);
    for i=linspace(0,1,speed)
        q_cur = quaternion.slerp(q1,q2,i,eps);
        [z,y,x]=quaternion.quat2angle(q_cur);
        l2 = spherical_harmonics.rotate(l1,...
            'xrotate',x,'yrotate',y,'zrotate',z);
        im2Vr = max(im1r(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vg = max(im1g(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vb = max(im1b(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2(:,:,1) = reshape(im2Vr,[size(im1,1) size(im1,2)] );
        im2(:,:,2) = reshape(im2Vg,[size(im1,1) size(im1,2)] );
        im2(:,:,3) = reshape(im2Vb,[size(im1,1) size(im1,2)] );
        imshow(im2);
        figure(f)
        drawnow
        frame = getframe();
        writeVideo(writerObj,frame);
        frames{end+1} = frame.cdata ;
    end
    q1 = q2;
    q2 = quaternion.angle2quat(0,-pi/4*span,pi/4*span);
    for i=linspace(0,1,speed*2)
        q_cur = quaternion.slerp(q1,q2,i,eps);
        [z,y,x]=quaternion.quat2angle(q_cur);
        l2 = spherical_harmonics.rotate(l1,...
            'xrotate',x,'yrotate',y,'zrotate',z);
        im2Vr = max(im1r(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vg = max(im1g(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vb = max(im1b(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2(:,:,1) = reshape(im2Vr,[size(im1,1) size(im1,2)] );
        im2(:,:,2) = reshape(im2Vg,[size(im1,1) size(im1,2)] );
        im2(:,:,3) = reshape(im2Vb,[size(im1,1) size(im1,2)] );
        imshow(im2);
        figure(f)
        drawnow
        frame = getframe();
        writeVideo(writerObj,frame);
        frames{end+1} = frame.cdata ;
    end
    q1 = q2;
    q2 = quaternion.angle2quat(0,-pi/4*span,-pi/4*span);
    for i=linspace(0,1,speed*2)
        q_cur = quaternion.slerp(q1,q2,i,eps);
        [z,y,x]=quaternion.quat2angle(q_cur);
        l2 = spherical_harmonics.rotate(l1,...
            'xrotate',x,'yrotate',y,'zrotate',z);
        im2Vr = max(im1r(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vg = max(im1g(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vb = max(im1b(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2(:,:,1) = reshape(im2Vr,[size(im1,1) size(im1,2)] );
        im2(:,:,2) = reshape(im2Vg,[size(im1,1) size(im1,2)] );
        im2(:,:,3) = reshape(im2Vb,[size(im1,1) size(im1,2)] );
        imshow(im2);
        figure(f)
        drawnow
        frame = getframe();
        writeVideo(writerObj,frame);
        frames{end+1} = frame.cdata ;
    end
    q1 = q2;
    q2 = quaternion.angle2quat(0,pi/4*span,-pi/4*span);
    for i=linspace(0,1,2*speed)
        q_cur = quaternion.slerp(q1,q2,i,eps);
        [z,y,x]=quaternion.quat2angle(q_cur);
        l2 = spherical_harmonics.rotate(l1,...
            'xrotate',x,'yrotate',y,'zrotate',z);
        im2Vr = max(im1r(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vg = max(im1g(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vb = max(im1b(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2(:,:,1) = reshape(im2Vr,[size(im1,1) size(im1,2)] );
        im2(:,:,2) = reshape(im2Vg,[size(im1,1) size(im1,2)] );
        im2(:,:,3) = reshape(im2Vb,[size(im1,1) size(im1,2)] );
        imshow(im2);
        figure(f)
        drawnow
        frame = getframe();
        writeVideo(writerObj,frame);
        frames{end+1} = frame.cdata ;
    end
    q1 = q2;
    q2 = quaternion.angle2quat(0,pi/4*span,0);
    for i=linspace(0,1,speed)
        q_cur = quaternion.slerp(q1,q2,i,eps);
        [z,y,x]=quaternion.quat2angle(q_cur);
        l2 = spherical_harmonics.rotate(l1,...
            'xrotate',x,'yrotate',y,'zrotate',z);
        im2Vr = max(im1r(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vg = max(im1g(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2Vb = max(im1b(:) .* max((l2'*Y)',0)./(l1'*Y)',0);
        im2(:,:,1) = reshape(im2Vr,[size(im1,1) size(im1,2)] );
        im2(:,:,2) = reshape(im2Vg,[size(im1,1) size(im1,2)] );
        im2(:,:,3) = reshape(im2Vb,[size(im1,1) size(im1,2)] );
        imshow(im2);
        figure(f)
        drawnow
        frame = getframe();
        writeVideo(writerObj,frame);
        frames{end+1} = frame.cdata ;
    end
    close(writerObj);
    vid_helper.make_gif(frames,'relighting.gif');
end

