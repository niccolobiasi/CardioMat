function lats=computeLats(filename,dt,thr,doplot)
% lats=computeLats(filename) returns a cell array whose each element is a
% numeric array containing the activation times for each voxel computed
% from the membrane potential.
%
% lats=computeLats(filename,dt) allows to specify the time resolution
% (default=1 ms).
%
% lats=computeLats(filename,dt,thr) allows to specify activation threshold
% (default=-20 mV).
%
% lats=computeLats(filename,dt,thr,doplot) allows to specify if a plot of
% the activation time should be generated. In this case a mat file
% containing the geometrical features should be in the current folder with
% the same name of the membrane potential file.


if nargin<2
    dt=1;
elseif isempty(dt)
    dt=1;
end


if nargin<3
    thr=-20;
elseif isempty(thr)
    thr=-20;
end

if nargin<4
    doplot=true;
elseif isempty(doplot)
    doplot=true;
end

load(filename,"ind_in");
sizeVm=nnz(ind_in);
fid=fopen(filename,'r');
Vm=fread(fid,sizeVm,'single');
prevVm=100*ones(sizeVm,1);
lats=cell(size(Vm));
i=1;
while ~isempty(Vm)
    ind_up=find(Vm>=thr & prevVm<thr);
    for j=1:length(ind_up)
        lats{ind_up(j)}=[lats{ind_up(j)} i*dt];
    end
    disp(['time =  ' num2str(i*dt)])
    prevVm=Vm;
    Vm=fread(fid,sizeVm,'single');
    i=i+1;
end
fclose(fid);
if doplot
    lat_first=extractLat(lats,'first');
    load(filename,"VoxelMat","plot_fib");
    latplot=NaN(size(VoxelMat));
    latplot(ind_in)=lat_first;
    [FV,extInd]=computeSurface(VoxelMat,1);
    figure;
    PlotVoxel(FV,latplot,extInd);
    hold on
    [FV1,extInd1]=computeSurface(plot_fib,1);
    PlotVoxel(FV1,plot_fib,extInd1,'k');
end
end