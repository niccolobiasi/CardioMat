function h=plotFrame(filename_bin,filename_geom, cmap)
%plotFrame(filename_bin) read the binary file with name filename_bin and
%allows interactive visualzation of saved data. plotFrame search for a mat
%file in the same path and with the same name containing simulation
%metadata (geometry adn fibrosis). Otherwise the mat file filename can be
%passed as second argument to the plotFrame function.

fid=fopen(filename_bin,'r');

if nargin<2
    filename_geom=filename_bin;
end
if nargin<3
    cmap=parula;
end
load(filename_geom,"VoxelMat","ind_in",'plot_fib','res');
sizeVm=nnz(ind_in);
Vm_plot=nan(size(VoxelMat));
fseek(fid,0,'eof');
Nframe=ftell(fid)/4/sizeVm;
fseek(fid,0,'bof');

[FV,myExtInd]=computeSurface(VoxelMat,res);

[FV1,myExtInd1]=computeSurface(plot_fib,res);


Vm=fread(fid,sizeVm,'single');
Vm_plot(ind_in)=Vm;


figure;
clf
h=PlotVoxel(FV,Vm_plot, myExtInd, cmap,1,0);
hold on
PlotVoxel(FV1,plot_fib,myExtInd1,'k',1,0);
colorbar
clim([-85 50])
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
frame=1;
disp(['Number of saved frame: ' num2str(Nframe)])
while true
    str=input('Type F to advance a frame, S to select a frame, Q to exit \n','s');
    if str=='Q'
        break;
    elseif str=='F'
        if frame==Nframe
            disp('Maximum frame achieved')
            continue;
        end
        Vm=fread(fid,sizeVm,'single');
        frame=frame+1;
    elseif str=='S'
        frame=input('Select frame:');
        if frame>Nframe || frame<1
            disp('Selcted frame does not exists')
            continue;
        end
        fseek(fid,4*sizeVm*(frame-1),'bof');
        Vm=fread(fid,sizeVm,'single');
    end

    Vm_plot(ind_in)=Vm;
    h.CData=Vm_plot(myExtInd);
    title(['frame=' num2str(frame)])

end
fclose(fid);
end

