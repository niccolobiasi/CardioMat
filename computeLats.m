function lats=computeLats(filename,dt,thr,start_time, end_time, doplot)
% lats=computeLats(filename) returns a cell array whose each element is a
% numeric array containing the activation times for each voxel computed
% from the membrane potential.
%
% lats=computeLats(filename,dt) allows to specify the time resolution
% (default=1 ms).
%
% lats=computeLats(filename,dt,thr,start_time, end_time) allows to specify
% activation threshold (default=-20 mV), and start and end time for
% scanning (default= 0 Inf)
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
    start_time=0;
end

if nargin<5
    end_time=Inf;
end


if nargin<6
    doplot=false;
elseif isempty(doplot)
    doplot=false;
end

load([filename '.mat'],"ind_in");
sizeVm=nnz(ind_in);
fid=fopen(filename,'r');
start_frame=round(start_time)/dt+1;
fseek(fid,4*sizeVm*(start_frame-1),'bof');
Vm=fread(fid,sizeVm,'single');
prevVm=100*ones(sizeVm,1);
lats=cell(size(Vm));
i=1;
while ~isempty(Vm) && i*dt+start_time<=end_time
    ind_up=find(Vm>=thr & prevVm<thr);
    for j=1:length(ind_up)
        lats{ind_up(j)}=[lats{ind_up(j)}; i*dt];
    end
    disp(['time =  ' num2str(i*dt)])
    prevVm=Vm;
    Vm=fread(fid,sizeVm,'single');
    i=i+1;
end
fclose(fid);
if doplot
    lat_first=extractLat(lats,'first');
    load([filename '.mat'],"VoxelMat","plot_fib");
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