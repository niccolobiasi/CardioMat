function runSimulation(VoxelMat,options)
% runSimulation(VoxelMat,options) runs a monodomain simulation on the
% voxelized geometry VoxelMat by using the smoothed boundary method.
% options is an optional structure input to define optional parameters.
% Possible fields of option structure are:
%
% - res: define the resolution of voxelized geometry (default=0.25 mm)
%
% - domain: can be a matrix of the same size of VoxelMat or a vector of
% size ninX1 with nin=nnz(VoxelMat). domain parameters define a label for
% each voxel in the geometry (default =1).
%
% - model: a string array defining the ionic model. The default model is
% Biasi_Tognetti. Other possible choice are:
%   - TNNP_2004
%   - CRN_1998
%   - TP_2006
%   - Custom model (see MyModel.m for a template function).
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
%       - CV (optional): conduction velocity in the Purkinje network
%       (default=340 cm/s)
%       - stimTimes (optional): time at which His bundle is stimulated
%       (default=0)
%       - current (optional): current delivered at the PMJ (default=300
%       V/s)
%       - duration (optional): duration of stimulation at the PMJ
%       (default=2 ms)
%       - ant_delay (optional): anterograde conduction delay at PMJ
%       (default=8 ms)
%       - retro_delay (optional): retrograde conduction delay at PMJ
%       (default=3 ms)
%       - ERP (optional): effective refractory period at PMJ(default=400ms)
%       - fib_thr (optional): thrshold value on the normalized pixel
%       intensity to disable PMJ in the border zone (default=0.4)
% duration, ant_delay, retro_delay, and ERP can be defined as scalar
% quantity equal for each PMJ or as a vector of length equal to terminal
% nodes in the purkinje tree (use "term_mask=not(ismember(pkn_elem(:,2),pkn_elem(:,1)));
% term_ind=pkn_elem(term_mask,2); " to find voxel coordinate of PMJs).
%
% - ExtStim: is an array of structures to define external stimulations with the
% following subfields (each array element correspond to a stimulation):
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
%binary file. The savefile structure must have the following subfields:
%       - filename: name of the file where the memebrane potential will be
%       saved
%       - dt: time interval for saving.
%
%- phase_xi: xi parameter used for phase field computation (see
%computePhase, default: 0.25 mm).
%
%- phase_dt: timestep used for phase field computation(default: 0.02 ms).
%
%- phase_thr: defines a phase field threshold below which voxels are not
%considered for computation of membrane potential (default: 1e-4).
%
%- save_state: structure to specify savings of the state of the
%simulations. The savestate structure may have two subfields:
%   - times: defines a vector of times at which saving the state
%     variables (default: [], no savings).
%   - filename: prefix name of mat files (default: 'sim').
%
%- load_state: filename of a mat file containing the state variables of the
%ionic model to use as initial condition.



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
    Stim(length(options.ExtStim))=struct();
    for i=1:length(options.ExtStim)
        Stim(i).ExtTimes=options.ExtStim(i).times;
        if isfield(options.ExtStim(i),'location')
            Stim(i).ExtLoc=options.ExtStim(i).location;
            Stim(i).toSelect=0;
        else
            Stim(i).toSelect=1;
        end
        Stim(i).Iext=options.ExtStim(i).current;
        Stim(i).ExtDuration=options.ExtStim(i).duration;
    end
else
    if isfield(options,'pkj') || isfield(options,'load_state')
        Stim(1).ExtTimes=[];
        Stim(1).ExtLoc=0;
        Stim(1).Iext=0;
        Stim(1).toSelect=0;
        Stim(1).ExtDuration=0;
    else
        Stim(1).ExtTimes=0;
        Stim(1).toSelect=1;
        Stim(1).Iext=50;
        Stim(1).ExtDuration=1;
    end
end

if isfield(options,'dt_visual')
    visual=round(options.dt_visual)/dt;
else
    visual=round(5/dt);
end

if isfield(options,'savefile')
    filename=options.savefile.filename;
    dt_save=options.savefile.dt;
    i_save=round(dt_save/dt);
    toSave=1;
else
    toSave=0;
end

if isfield(options,'phase_xi')
    xi=options.phase_xi/10;
else
    xi=0.025;
end

if isfield(options,'phase_dt')
    phase_dt=options.phase_dt;
else
    phase_dt=0.02;
end

if isfield(options,'phase_thr')
    thr_b=options.phase_thr;
else
    thr_b=1e-4;
end

if isfield(options,'save_state')
    if isfield(options.save_state,'times')
        savestate_times=options.save_state.times;
    else
        savestate_times=[];
    end
    if isfield(options.save_state,'filename')
        filename_state=options.save_state.filename;
    elseif isfield(options.savefile,'filename')
        filename_state=options.savefile.filename;
    else
        filename_state='sim';
    end
else
    savestate_times=[];
end

if isfield(options,'model')
    model=options.model;
else
    model='Biasi_Tognetti';
end

%% Define simulation length
t=gpuArray(0:dt:T);
Nt=length(t);

%% Phase field and mask computation

phi=computePhase(VoxelMat,res/10, phase_dt,xi,true);
mask_phi=phi>thr_b;


%% extrapolate fibers and domain and pixel intensity
[~,extInd]=computeSurface(VoxelMat,res);
fields=extrapField([f.x,f.y,f.z,domain,scalar,plot_fib(VoxelMat)],VoxelMat,extInd,mask_phi);
fx=fields(:,1);
fy=fields(:,2);
fz=fields(:,3);
fib=gpuArray(fields(:,6));
plot_fib2=false(size(VoxelMat));
plot_fib2(mask_phi)=logical(fib(mask_phi(:)));
ind_in= mask_phi & not(plot_fib2);
nin=nnz(ind_in);
domain=gpuArray(fields(ind_in(:),4));
scalar=gpuArray(fields(ind_in(:),5));

clear fields
%% system matrix
As=gpuArray(heartMatrixFib(phi,mask_phi,plot_fib2,fx,fy,fz,Dpar,anisot_ratio,res));
%% elab purkinje system

if isfield(options,'pkj')
    pkn_nodes=options.pkj.nodes;
    pkn_elem=options.pkj.elem;
    term_mask=not(ismember(pkn_elem(:,2),pkn_elem(:,1)));
    term_ind=pkn_elem(term_mask,2);
    term_ind=[1;sort(term_ind)];
    distM=dist_pkj(pkn_nodes,pkn_elem);
    if isfield(options.pkj,'CV')
        pkj_vel=options.pkj.CV;
    else
        pkj_vel=340;
    end
    distM=distM(term_ind,term_ind)*res/pkj_vel*100;
    term_ind=term_ind(2:end);
    term_nodes=pkn_nodes(term_ind,:);
    term_map=round(sub2ind(size(VoxelMat),term_nodes(:,2),term_nodes(:,1),term_nodes(:,3)));
    tmp=zeros(size(VoxelMat));
    tmp(term_map)=2:length(term_nodes)+1;
    tmp=tmp(ind_in);
    if isfield(options.pkj,'fib_thr')
        fib_thr=options.pkj.fib_thr;
    else
        fib_thr=0.4;
    end
    tmp(scalar>fib_thr)=0;
    term_map2=find(tmp);
    good_term=tmp(term_map2);
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
    if isfield(options.pkj,'ant_delay')
        ant_delay=options.pkj.ant_delay;
    else
        ant_delay=8;
    end
    if isfield(options.pkj,'retro_delay')
        retro_delay=options.pkj.retro_delay;
    else
        retro_delay=3;
    end
    if isfield(options.pkj,'ERP')
        ERP=options.pkj.ERP;
    else
        ERP=400;
    end
else
    pkn_stim=[];
end
pkn_stim_n=round(pkn_stim/dt)+1;

%% Set up external stimulation
for i=1:length(Stim)
    if Stim(i).toSelect
        disp(['Select stimulation region ', num2str(i), ' by mouse click']);
        Stim(i).ExtLoc=selectVolume(VoxelMat,res,5,false, plot_fib);
    end
    I_mask=zeros(size(VoxelMat));
    I_mask(Stim(i).ExtLoc & VoxelMat)=Stim(i).Iext;
    Stim(i).I_mask=gpuArray(I_mask(ind_in));

    n_stim=length(Stim(i).ExtTimes);

    NT=round(Stim(i).ExtDuration/dt)+1;

    Stim(i).ext_stim=false(Nt,1);
    for j=1:n_stim
        ind=round(Stim(i).ExtTimes(j)/dt)+1;
        Stim(i).ext_stim(ind:ind+NT)=1;
    end
end
%% Initialization
[Vm0, state_model0]=feval(model);
nVar=length(state_model0);
if isfield(options,'load_state')
    saved_state=load(options.load_state);
    Vm=saved_state.Vm;
    state_model=saved_state.state_model;
    if length(state_model)~=nVar
        error('number of loaded state variables do not correspond to the selected model');
    end
else
    Vm=gpuArray(Vm0*ones(nin,1));
    state_model=cell(nVar,1);
    for ii=1:nVar
        state_model{ii}=state_model0{ii}*ones(nin,1);
    end
end

Ie=gpuArray(zeros(nin,1));

if ~isinf(visual)
    Vm_plot=NaN(size(VoxelMat));
    Vm_plot(ind_in)=gather(Vm);
    [FV1,extInd1]=computeSurface(plot_fib,res);
    hf=figure;
    hf.MenuBar='figure';
    [FV,extInd]=computeSurface(VoxelMat & ~plot_fib,res);
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
    save(filename,'ind_in','plot_fib','VoxelMat','res','mask_phi','dt_save');
end
%% initializa PKJ
if  exist('term_map2','var')
    state=gpuArray(zeros(size(good_term)));
    t_last=gpuArray(zeros(size(good_term)));
    Tact=Inf(size(good_term))';
end

%% RUN SIMULATION

for i=1:Nt
    dv2dt=As*Vm;
    for j=1:length(Stim)
        if Stim(j).ext_stim(i)
            dv2dt=dv2dt+Stim(j).I_mask;
        end
    end

    [Vm,state_model]=feval(model,Vm,state_model,Ie,dv2dt,domain,dt,scalar);

    if toSave
        if mod(i,i_save)==0
            fwrite(fid,gather(Vm),'single');
        end
    end

    if ismember(dt*i,savestate_times)
        save([filename_state '_state_' num2str(dt*i)], 'Vm','state_model');
    end

    if exist('term_map2','var')
        active=dv2dt(term_map2)>0 & Vm(term_map2)>-20;
        if ismember(i,pkn_stim_n)
            Tact=min(Tact,distM(1,good_term));
        end
        Tact=Tact-dt;
        [s_v,s_pkj,state,t_last]=PKJnodeFun(Tact'<0,active,t(i),state,t_last,ant_delay,retro_delay,pkj_duration,0,ERP);
        Tact(Tact<0)=Inf;
        Tact=min([Tact;distM(good_term(s_pkj),good_term)],[],1);
        Ie(:)=0;
        Ie(term_map2(s_v))=pkj_current;
    end
    if mod(i,visual)==0
        Vm_plot(ind_in)=gather(Vm);
        h.CData=Vm_plot(extInd);
        title(num2str(dt*i))
        drawnow
    end

    if exist('term_map2','var') 
         check_pkj=nnz(~isinf(Tact))==0 && nnz(state<3 & state>0)==0 && nnz(pkn_stim(i:end))==0;
    else 
        check_pkj=1;
    end

    if nnz(Vm>=-40)<5 && check_pkj
        stop=1;
        for j=1:length(Stim)
            if nnz(Stim(j).ext_stim(i:end))>0
                stop=0;
                break;
            end
        end
        if stop==1
            break;
        end
    end
end

if toSave
    fclose(fid);
end


