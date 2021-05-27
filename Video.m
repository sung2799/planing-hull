% Video
v = VideoWriter('trajectory.avi','Motion JPEG AVI');
v = VideoWriter('trajectory.avi');
v.FrameRate=15;
open(v);

for k=1:length(t)
    drawroll(yNL(k,:));
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(.01);
end
close(v);
