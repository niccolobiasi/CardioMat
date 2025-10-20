function createVideo_wh(filename_save,filename_bin,filename_mat,video_res,plane,show_fib)
% createVideo_wh is the equivalent of createVideo for whole heart 
% simulations. It reads the binary file filename_bin and creates a mp4
% video of the membrane potential map. The function assumes there is a mat
%  file containing simulation metadata in the same path and with the same
%  name of the binary file. If not, the user can specify the mat file.
% Additionally, the user can specify the video resolution video_res and a
% cutplane defined as a logical matrix of the same size of the
% computational box(see DefineCP.m). If show_fib is set to 0 the fibrosis
% is not shown (default=1).

if nargin<3
    filename_mat=[filename_bin '.mat'];
end
if isempty(filename_mat)
    filename_mat=[filename_bin '.mat'];
end


fid=fopen(filename_bin,'r');
load(filename_mat,"ventr","atr","ind_in_atr",'ind_in_ventr','plot_fib','res','dt_save');
VoxelMat=atr | ventr;
sizeVm_atr=nnz(ind_in_atr);
sizeVm_ventr=nnz(ind_in_ventr);
Vm_plot=nan(size(VoxelMat));

if nargin<4
    video_res=150;
end

if isempty(video_res)
    video_res=150;
end

if nargin<5
    plane=true(size(VoxelMat));
end

if isempty(plane)
    plane=true(size(VoxelMat));
end

if nargin<6
    show_fib=1;
end

[FV,extInd]=computeSurface(VoxelMat & ~plot_fib & plane,res);
if show_fib
    [FV1,myExtInd1]=computeSurface(plot_fib & plane,res);
end

Vm=fread(fid,sizeVm_atr,'single');
Vm_plot(ind_in_atr)=Vm;
Vm=fread(fid,sizeVm_ventr,'single');
Vm_plot(ind_in_ventr)=Vm;


hf=figure;
hf.MenuBar='figure';
clf
h=PlotVoxel(FV,Vm_plot, extInd, parula,1,0);
if show_fib
    hold on
    PlotVoxel(FV1,plot_fib,myExtInd1,'k',1,0);
end
colorbar
clim([-85 50])
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
view(3)
ttl=title('Time=0 ms');
drawnow

while true
    str=input('Orient geometry and light for the video, once done type Q \n','s');
    if str=='Q'
        break;
    end
end

video=VideoWriter(filename_save,'MPEG-4');
open(video);
RGBImage=print('-RGBImage',['-r' num2str(video_res)]);
frame=im2frame(RGBImage);
writeVideo(video,frame);


i=2;
Vm_atr=fread(fid,sizeVm_atr,'single');
Vm_ventr=fread(fid,sizeVm_ventr,'single');
while ~isempty(Vm_atr) && ~isempty(Vm_ventr)

    Vm_plot(ind_in_atr)=Vm_atr;
    Vm_plot(ind_in_ventr)=Vm_ventr;
    h.CData=Vm_plot(extInd);
    ttl.String=['Time=' num2str((i-1)*dt_save)];
    drawnow

    RGBImage=print('-RGBImage',['-r' num2str(video_res)]);
    frame=im2frame(RGBImage);
    writeVideo(video,frame);

    i=i+1;
    Vm_atr=fread(fid,sizeVm_atr,'single');
    Vm_ventr=fread(fid,sizeVm_ventr,'single');
end
fclose(fid);
close(video);
end

