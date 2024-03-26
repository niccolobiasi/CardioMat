function runSimulation(VoxelMat,options)
% runSimulation(VoxelMat,options) runs a monodomain simulation on the
% voxelized geometry VoxelMat by using the smoothed boundary method.
% options is an optional structure input to define optional parameters.
% Possible field of options structures are:
%
% - res: define the resolution of voxelized geometry (default=0.25 mm)
%
% - domain: can be a matrix of the same size of VoxelMat or a vector of
% size ninX1 with nin=nnz(VoxelMat). domain parameters define a label for
% each voxel in the geometry (default =1).
%
% - fib: is a matrix of the same size of VoxelMat indicating fibrotic\scar
% voxels (default=0);
%
% - scalar: can be a matrix of the same size of VoxelMat or a vector of
% size ninX1. scalar defines the normalized pixel intensity to be used in
% the ionic model for electrophysiological adjustments (default=0).
%
% - fibers: is a structure with three fields ('x','y','z') containing the
% three components of fiber direction (default=1,0,0). Each field must be a
% vector of size ninX1.
%
% - Dpar: defines the diffusion in the longitudinal fiber direction
% (default= 1.2e-3)
%
% - anisot_ratio: defines the anisotropy ratio (default=4)
%
% - dt: time step used for forward Euler method (default=0.02 ms)
%
% - Tsimu: total simulation time (default= 200 ms)
%
% - pkj: define the purkinje network. It is a structure with the following
% subfields:
%       - nodes (mandatory): nodes composing Purkinje network
%       - elem (mandatory): nodes composing Purkinje network
%       - times (mandatory): activation times for each node of PKJ network
%       - stimTimes (optional): time at which His bundle is stimulated
%       (default=0)
%       - current (optional): current delivered at the PMJ (default=300
%       V/s)
%       - duration (optional): duration of stimulation at the PMJ
%       (default=2 ms)
%
% - ExtStim: is a structure to define an external stimulation with the
% following subfields:
%       - times (mandatory): a vector containing the stimulation times
%       - duration (mandatory): duration of the stimulations
%       - current (mandatory): stimulation current
%       - location (optional): a logical matrix of the same size of VoxelMat
%       defining the region to stimulate. If location is not provided the
%       user is asked to select the regions to stimulate by mouse click.
%  In the case no both ExtStim and pkj fields are not passed to the
%  function, ExtStim will be automatically set up with a single stimulation
%  at the begining of the simulation with current equal to 50 and duration
%  equal to 1. The location will be selected by the user.
%
% - dt_visual: defines time interval for visualization. If dt_visual=Inf
% then no visualization is provided (default= 5 ms).
%
%- savefile: structure to specify if save the membrane potential to a
%binary file. the savefile structure must have the following subfields:
%       - filename: name of the file where the memebrane potential will be
%       saved
%       - dt: time interval for saving.




if nargin<2
    options=[];
end

nin=nnz(VoxelMat);

if isfield(options,'domain')
    domain=options.domain;
    if numel(domain)==numel(VoxelMat)
        domain=domain(VoxelMat(:));
    end
else
    domain=ones(nin,1);
end

if isfield(options,'fib')
    plot_fib=options.fib;
else
    plot_fib=false(size(VoxelMat));
end

if isfield(options,'scalar')
    scalar=options.scalar;
    if numel(scalar)==numel(VoxelMat)
        scalar=scalar(VoxelMat(:));
    end
else
    scalar=zeros(nin,1);
end

if isfield(options,'res')
    res=options.res;
else
    res=0.25;
end

if isfield(options,'fibers')
    f=options.fibers;
else
    f.x=ones(nin,1);
    f.y=zeros(nin,1);
    f.z=zeros(nin,1);
end

if isfield(options,'Dpar')
    Dpar=options.Dpar;
else
    Dpar=1.2e-3;
end


if isfield(options,'anisot_ratio')
    anisot_ratio=options.anisot_ratio;
else
    anisot_ratio=4;
end

if isfield(options,'dt')
    dt=options.dt;
else
    dt=0.02;
end


if isfield(options,'Tsimu')
    T=options.Tsimu;
else
    T=200;
end

if isfield(options, 'ExtStim')
    ExtTimes=options.ExtStim.times;
    if isfield(options.ExtStim,'locations')
        ExtLoc=options.ExtStim.location;
        toSelect=0;
    else
        toSelect=1;
    end
    Iext=options.ExtStim.current;
    ExtDuration=options.ExtStim.duration;
else
    if isfield(options,'pkj')
        ExtTimes=[];
        ExtLoc=0;
        Iext=0;
        toSelect=0;
        ExtDuration=0;
    else
        ExtTimes=0;
        toSelect=1;
        Iext=50;
        ExtDuration=1;
    end
end

if isfield(options,'dt_visual')
    visual=round(options.dt_visual)/dt;
else
    visual=round(5/dt);
end

if isfield(options,'savefile')
    filename=options.savefile.filename;
    i_save=round(options.savefile.dt/dt);
    toSave=1;
else
    toSave=0;
end

%% Define simulation length
t=gpuArray(0:dt:T);
Nt=length(t);

%% Phase field and mask computation

phi=computePhase(VoxelMat,res/10, 0.02,0.025,true);
thr_b=1e-4;
mask_phi=phi>thr_b;
ind_in= mask_phi & not(plot_fib);
nin=nnz(ind_in);

%% extrapolate fibers and domain and pixel intensity
[FV,extInd]=computeSurface(VoxelMat,res);
fields=extrapField([f.x,f.y,f.z,domain,scalar],VoxelMat,extInd,mask_phi);
fx=fields(:,1);
fy=fields(:,2);
fz=fields(:,3);
domain=gpuArray(fields(ind_in(:),4));
scalar=gpuArray(fields(ind_in(:),5));
clear fields
%% system matrix
As=gpuArray(heartMatrixFib(phi,mask_phi,plot_fib,fx,fy,fz,Dpar,anisot_ratio,res));
%% elab purkinje system

if isfield(options,'pkj')
    pkn_nodes=options.pkj.nodes;
    pkn_elem=options.pkj.elem;
    pkn_times=options.pkj.times;
    term_mask=not(ismember(pkn_elem(:,2),pkn_elem(:,1)));
    term_ind=pkn_elem(term_mask,2);
    term_nodes=pkn_nodes(term_ind,:);
    term_times=pkn_times(term_ind);
    term_map=sub2ind(size(VoxelMat),term_nodes(:,2),term_nodes(:,1),term_nodes(:,3));
    tmp=zeros(size(VoxelMat));
    tmp(term_map)=1:length(term_nodes);
    tmp=tmp(ind_in);
    term_map2=find(tmp);
    good_term=tmp(term_map2);
    act_times=term_times(good_term);
    Tact=Inf(size(good_term));
    clear tmp
    if isfield(options.pkj,'stimTimes')
        pkn_stim=options.pkj.stimTimes;
    else
        pkn_stim=0;
    end
    if isfield(options.pkj,'current')
        pkj_current=options.pkj.current;
    else
        pkj_current=300;
    end
    if isfield(options.pkj,'duration')
        pkj_duration=options.pkj.duration;
    else
        pkj_duration=2;
    end
else
    pkn_stim=[];
end
pkn_stim_n=round(pkn_stim/dt)+1;

%% Set up external stimulation
if toSelect
    disp('Select stimulation regions by mouse click')
    ExtLoc=selectVolume(VoxelMat,res,5,false, plot_fib);
end
I_mask=zeros(size(VoxelMat));
I_mask(ExtLoc & VoxelMat)=Iext;
I_mask=gpuArray(I_mask(ind_in));

n_stim=length(ExtTimes);
ind=0;
NT=round(ExtDuration/dt)+1;

ext_stim=false(Nt,1);
for j=1:n_stim
    ind=ind+round(ExtTimes(j)/dt)+1;
    ext_stim(ind:ind+NT)=1;
end
%% Initialization

Vm=gpuArray(-85*ones(nin,1));
u=gpuArray(zeros(nin,1));
w=gpuArray(zeros(nin,1));

Ie=gpuArray(zeros(nin,1));

if ~isinf(visual)
    Vm_plot=NaN(size(VoxelMat));
    Vm_plot(ind_in)=gather(Vm);
    [FV1,extInd1]=computeSurface(plot_fib,res);
    figure;
    h=PlotVoxel(FV,Vm_plot, extInd, parula,1,0);
    hold on
    PlotVoxel(FV1,plot_fib,extInd1,'k',1,0);
    colorbar
    clim([-85 50]);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    camlight
end
%% create file for saving

if toSave
    fid=fopen(filename,'w+');
    save(filename,'ind_in','plot_fib','VoxelMat','res');
end
%% RUN SIMULATION

for i=1:Nt-1
    dv2dt=As*Vm;
    if ext_stim(i)
        dv2dt=dv2dt+I_mask;
    end
    if ~isempty(pkn_stim)
        if nnz(i==pkn_stim_n)
            Tact=min(Tact,act_times);
        end
        Tact=Tact-dt;
        Ie(:)=0;
        Ie(term_map2(Tact<0))=pkj_current;
        Tact(Tact<-pkj_duration)=Inf;
    end
    [Vm,u,w]=single_step(Vm,u,w,Ie,dv2dt,domain,dt,scalar);
    if toSave
        if mod(i,i_save)==0
            fwrite(fid,gather(Vm),'single');
        end
    end

    if mod(i,visual)==0
        Vm_plot(ind_in)=gather(Vm);
        h.CData=Vm_plot(extInd);
        title(num2str(dt*i))
        drawnow
    end
end

if toSave
    fclose(fid);
end


