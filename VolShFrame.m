function vol_handle=VolShFrame(filename_bin,filename_geom,cmap)
%VolShowFrame(filename_bin) read the binary file with name filename_bin and
%allows interactive visualzation of saved data by using MATLAB volume. 
% VolShowFrame search for a mat file in the same path and with the same
% name containing simulation metadata. Otherwise the mat file filename can be
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
fseek(fid,0,'eof');
Nframe=ftell(fid)/4/sizeVm;
fseek(fid,0,'bof');

ind_nofib=VoxelMat & ~plot_fib;

[~,extInd]=computeSurface(ind_nofib,res);


Vm=fread(fid,sizeVm,'single');


idx=reshape(extrapField((1:length(Vm))',ind_in,extInd,[],10),size(VoxelMat));
idx_toadj=~isnan(idx) & ~ind_in;
idx=idx(idx_toadj);
vol_col=zeros(size(VoxelMat));
vol_col(ind_in)=Vm;
vol_col(idx_toadj)=Vm(idx);

vol_handle=volshow(vol_col);
vol_handle.AlphaData=ind_nofib;
vol_handle.Colormap=cmap;
viewer=vol_handle.Parent;
viewer.BackgroundColor=[1 1 1];
viewer.BackgroundGradient='off';
vol_fib=volshow(plot_fib,Parent=viewer);
vol_fib.Colormap=[0 0 0];
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

    vol_handle.Data(ind_in)=Vm;
    vol_handle.Data(idx_toadj)=Vm(idx);
    disp(['frame=' num2str(frame)])

end

end

