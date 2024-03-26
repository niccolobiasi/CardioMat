function createVideo(filename_save,filename_bin,filename_mat,dt)
%createVideo(filename_save,filename_bin) reads the binary file filename_bin
%and creates a mp4 video of the membrane potential map. The function
%assumes there is a mat file containing simulation metadata in the same
%path and with the same name of the binary file. If not, the user can specify
% the mat file. Additionally, the user can specify the time interval 
% between frames dt


if nargin<3
    filename_mat=filename_bin;
end
if isempty(filename_mat)
    filename_mat=filename_bin;
end
if nargin<4
    dt=1;
end
fid=fopen(filename_bin,'r');
load(filename_mat,'VoxelMat',"ind_in",'plot_fib','res');
sizeVm=nnz(ind_in);
Vm_plot=nan(size(VoxelMat));



[FV,myExtInd]=computeSurface(VoxelMat,res);

[FV1,myExtInd1]=computeSurface(plot_fib,res);


Vm=fread(fid,sizeVm,'single');
Vm_plot(ind_in)=Vm;


figure;
clf
h=PlotVoxel(FV,Vm_plot, myExtInd, parula,1,0);
hold on
PlotVoxel(FV1,plot_fib,myExtInd1,'k',1,0);
colorbar
clim([-85 50])
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
view(3)
title(num2str(dt))
drawnow

while true
    str=input('Orient geometry and light for the video, once done type Q \n','s');
    if str=='Q'
        break;
    end
end

video=VideoWriter(filename_save,'MPEG-4');
open(video);
RGBImage=print('-RGBImage','-r300');
frame=im2frame(RGBImage);
writeVideo(video,frame);


i=2;
Vm=fread(fid,sizeVm,'single');
while ~isempty(Vm)

    Vm_plot(ind_in)=Vm;
    h.CData=Vm_plot(myExtInd);
    title(num2str(dt*i))
    drawnow

    RGBImage=print('-RGBImage','-r300');
    frame=im2frame(RGBImage);
    writeVideo(video,frame);

    i=i+1;
    Vm=fread(fid,sizeVm,'single');

end
fclose(fid);
close(video);
end

