function make_vid( s,fname)
%MAKE_VID Summary of this function goes here
%   Detailed explanation goes here
set(gcf,'Color',[1 1 1]);
axis off

writerObj = VideoWriter(fname,'MPEG-4');
writerObj.FrameRate = 30;
open(writerObj);

view([0 -70])
rotate(s,[0 1 0],-85)
for i=1:120
    rotate(s,[0 1 0],4/3)
    writeVideo(writerObj,getframe);
end

for i=-70:-1.3:-90
    view([0 i])
    writeVideo(writerObj,getframe);
end
for i=-91:-1.3:-110
    view([0 i])
    camroll(180)
    writeVideo(writerObj,getframe);
end

for i=1:120
    rotate(s,[0 1 0],-4/3)
    writeVideo(writerObj,getframe);
end

for i=-109:1.3:-91
    view([0 i])
    camroll(180)
    writeVideo(writerObj,getframe);
end
for i=-90:1.3:-70
    view([0 i])
    writeVideo(writerObj,getframe);
end
close(writerObj);

end

