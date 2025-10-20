function [vol_handle_atr, vol_handle_ventr, viewer]=VolShFrame_wh(filename_bin,filename_geom,cmap)
% VolShowFrame_wh(filename_bin) is th equivalent of VolShowFrame for whole
% heart simulations. It reads the binary file with name filename_bin and
% allows interactive visualzation of saved data by using MATLAB volume. 
% VolShowFrame_wh search for a mat file in the same path and with the same
% name containing simulation metadata. Otherwise the mat file filename can 
% be passed as second argument.

fid=fopen(filename_bin,'r');

if nargin<2
    filename_geom=filename_bin;
end
if nargin<3
    cmap=parula;
end

load(filename_geom,"ventr","atr","ind_in_atr",'ind_in_ventr','plot_fib','res','dt_save');
VoxelMat=atr | ventr;
sizeVm_atr=nnz(ind_in_atr);
sizeVm_ventr=nnz(ind_in_ventr);
fseek(fid,0,'eof');
Nframe=ftell(fid)/4/(sizeVm_atr+sizeVm_ventr);
fseek(fid,0,'bof');

ind_nofib_atr=atr & ~plot_fib;
ind_nofib_ventr=ventr & ~plot_fib;

[~,extInd_atr]=computeSurface(ind_nofib_atr,res);
[~,extInd_ventr]=computeSurface(ind_nofib_ventr,res);


Vm_atr=fread(fid,sizeVm_atr,'single');
Vm_ventr=fread(fid,sizeVm_ventr,'single');

[Ny,Nx,~]=size(VoxelMat);

enlarge_atr=ind_nofib_atr;
enlarge_atr(extInd_atr+1)=1;
enlarge_atr(extInd_atr-1)=1;
enlarge_atr(extInd_atr+Ny)=1;
enlarge_atr(extInd_atr-Ny)=1;
enlarge_atr(extInd_atr+Ny*Nx)=1;
enlarge_atr(extInd_atr-Ny*Nx)=1;

enlarge_ventr=ind_nofib_ventr;
enlarge_ventr(extInd_ventr+1)=1;
enlarge_ventr(extInd_ventr-1)=1;
enlarge_ventr(extInd_ventr+Ny)=1;
enlarge_ventr(extInd_ventr-Ny)=1;
enlarge_ventr(extInd_ventr+Ny*Nx)=1;
enlarge_ventr(extInd_ventr-Ny*Nx)=1;

idx_atr=reshape(extrapField((1:length(Vm_atr))',ind_in_atr,extInd_atr,enlarge_atr),size(VoxelMat));
idx_toadj_atr=enlarge_atr & ~ind_in_atr;
idx_atr=idx_atr(idx_toadj_atr);

vol_col_atr=zeros(size(VoxelMat));
vol_col_atr(ind_in_atr)=Vm_atr;
vol_col_atr(idx_toadj_atr)=Vm_atr(idx_atr);


idx_ventr=reshape(extrapField((1:length(Vm_ventr))',ind_in_ventr,extInd_ventr,enlarge_ventr),size(VoxelMat));
idx_toadj_ventr=enlarge_ventr & ~ind_in_ventr;
idx_ventr=idx_ventr(idx_toadj_ventr);

vol_col_ventr=zeros(size(VoxelMat));
vol_col_ventr(ind_in_ventr)=Vm_ventr;
vol_col_ventr(idx_toadj_ventr)=Vm_ventr(idx_ventr);


vol_handle_atr=volshow(vol_col_atr);
vol_handle_atr.AlphaData=ind_nofib_atr;
vol_handle_atr.Colormap=cmap;
vol_handle_atr.DataLimits=[-90 50];
viewer=vol_handle_atr.Parent;
viewer.BackgroundColor=[1 1 1];
viewer.BackgroundGradient='off';
vol_fib=volshow(plot_fib,Parent=viewer);
vol_fib.Colormap=[0 0 0];
vol_fib.DataLimits=[0 1];
vol_handle_ventr=volshow(vol_col_ventr,Parent=viewer);
vol_handle_ventr.AlphaData=ind_nofib_ventr;
vol_handle_ventr.Colormap=cmap;
vol_handle_ventr.DataLimits=[-90 50];
frame=1;
disp(['Number of saved frame: ' num2str(Nframe)])
disp(['frame=' num2str(frame) ', t=' num2str((frame-1)*dt_save) 'ms'])
while true
    str=input('Type F to advance a frame, S to select a frame, Q to exit \n','s');
    if str=='Q'
        break;
    elseif (strcmp(str,'F') || strcmp(str,''))
        if frame==Nframe
            disp('Maximum frame achieved')
            continue;
        end
        Vm_atr=fread(fid,sizeVm_atr,'single');
        Vm_ventr=fread(fid,sizeVm_ventr,'single');
        frame=frame+1;
    elseif str=='S'
        frame=input('Select frame:');
        if frame>Nframe || frame<1
            disp('Selcted frame does not exists')
            continue;
        end
        fseek(fid,4*(sizeVm_atr+sizeVm_ventr)*(frame-1),'bof');
        Vm_atr=fread(fid,sizeVm_atr,'single');
        Vm_ventr=fread(fid,sizeVm_ventr,'single');
    end

    vol_handle_atr.Data(ind_in_atr)=Vm_atr;
    vol_handle_atr.Data(idx_toadj_atr)=Vm_atr(idx_atr);
    vol_handle_ventr.Data(ind_in_ventr)=Vm_ventr;
    vol_handle_ventr.Data(idx_toadj_ventr)=Vm_ventr(idx_ventr);
    disp(['frame=' num2str(frame) ', t=' num2str((frame-1)*dt_save) 'ms'])

end

end

