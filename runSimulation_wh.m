function runSimulation_wh(atr,ventr,options)
% runSimulation(atr,ventr,options) runs a monodomain simulation on the
% whole-heart geometry defined by atr and ventr matrices.
% atr and ventr are the atrial and ventricular voxelized geometries, 
% respectively. The numerical methods adopted are the same of runSimulation
% function. atr and ventr must be matrices of the same size.
% options is an optional structure input to define optional parameters.
% Possible fields of option structure are:
%
% - res: define the resolution of voxelized geometry (default=0.25 mm)
%
% - domain_atr: can be a matrix of the same size of atr or a vector of
% size nin_atrX1 with nin_atr=nnz(atr). domain_atr defines a label for
% each voxel in the atrial geometry (default =1).
%
% - domain_ventr: can be a matrix of the same size of ventr or a vector of
% size nin_ventrX1 with nin_ventr=nnz(ventr). domain_ventr defines a label 
% for each voxel in the ventricular geometry (default =1).
%
% - domain: a matrix of the same size of atr and ventr defining domains in
% both atrial and ventricular geometries (default=1). If this parameter is 
% passed domain_ventr and domain_atr are ignored.
%
% - model_atr: a string array defining the ionic model to be used in atrial
% geometry. The default model is Biasi_Tognetti. See runSimulation help for
% other possible choices.
%
% - model_ventr: a string array defining the ionic model to be used in 
% ventricular geometry. The default model is Biasi_Tognetti. See 
% runSimulation help for other possible choices.
%
% - fib_atr: can be a matrix of the same size of atr or a vector of
% size nin_atrX1 with nin_atr=nnz(atr). domain_atr defines fibrotic voxels
% in the atrial geometry (default =0).
%
% - fib_ventr: can be a matrix of the same size of ventr or a vector of
% size nin_ventrX1 with nin_ventr=nnz(ventr). domain_ventr defines fibrotic
% voxels in the ventricular geometry (default =0).
%
% - fib: is a matrix of the same size of atr and ventr indicating 
% fibrotic\scar voxels in both atria an ventricular geometries(default=0).
% If this parameter is passed fib_ventr and fib_atr are ignored.
%
% - scalar_atr: can be a matrix of the same size of atr or a vector of
% size nin_atrX1 with nin_atr=nnz(atr). domain_atr defines normalized pixel
% intensity to be used in the atrial ionic model for electrophysiological 
% adjustments (default=0).
%
% - scalar_ventr: can be a matrix of the same size of ventr or a vector of
% size nin_ventrX1 with nin_ventr=nnz(ventr). domain_ventr defines 
% normalized pixel intensity to be used in the ventricular ionic model for
% electrophysiological adjustments (default=0).
%
% - scalar: is a matrix of the same size of atr and ventr defining the 
% normalized pixel intensity to be used in both the atrial and ventricular
% ionic model for electrophysiological adjustments (default=0). If this 
% parameter is passed scalar_ventr and scalar_atr are ignored.
%
% - fibers_atr: is a structure with three fields ('x','y','z') containing 
% the three components of atrial fiber direction (default=1,0,0). Each 
% field must be a vector of size nin_atrX1.
%
% - fibers_ventr: is a structure with three fields ('x','y','z') containing 
% the three components of ventricular fiber direction (default=1,0,0). Each 
% field must be a vector of size nin_ventrX1.
%
% - Dpar_atr: defines the diffusion in the longitudinal fiber direction
% in the atrial geometry (default= 4.8e-3 cm^2/ms).
%
% - Dpar_ventr: defines the diffusion in the longitudinal fiber direction
% in the ventricular geometry (default= 1.2e-3 cm^2/ms).
%
% - anisot_ratio_atr: defines the anisotropy ratio in the atrial geometry
% (default=8).
%
% - anisot_ratio_ventr: defines the anisotropy ratio in the ventricular
% geometry (default=4).
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
%       - rv_delay: additional delay in the connection between His bundle
%       and right bundle branch (default=0ms)
%       - lv_delay: additional delay in the connection between His bundle
%       and left bundle branch (default=0ms)
%       - avn_anterograde_delay: anterograde conduction delay at the AVN
%       (default=80ms)
%       - avn_retrograde_delay: retrograde conduction delay at the AVN
%       (default=80ms)
%       - avn_stim_duration: duration of the atrial stimulation in
%       retrograde conduction (default=2ms);
%       - avn_current: atrial stimulation current in retrograde conduction
%       (default=200mV/ms)
%       - avn_erp: effective refrcatory period of the AVN (default=300ms)
%       - av_block: a string defining type of block at the level of the
%       AVN. Possible choices are "none"(default), "anterograde",
%       "retrograde", "both".
% duration, ant_delay, retro_delay, and ERP can be defined as scalar
% quantity equal for each PMJ or as a vector of length equal to terminal
% nodes in the purkinje tree (use "term_mask=not(ismember(pkn_elem(:,2),pkn_elem(:,1)));
% term_ind=pkn_elem(term_mask,2); " to find voxel coordinate of PMJs).
% The rv_delay and lv_delay parameters are applied by assuming the Purkinje
% network is generated with the createPurkinje_biv and insertAVN functions 
% (the number of node must be even with the first node and second nodes
% defining the AVN and His bundle, respectively; left and right purkinje 
% networks must be composed by the same number of nodes; the left purkinje
% nodes shoud appear first in the pkn_nodes matrix).
%
% - ExtStim: is an array of structures to define external stimulations with
% the following subfields (each array element correspond to a stimulation):
%       - times (mandatory): a vector containing the stimulation times
%       - duration (mandatory): duration of the stimulations
%       - current (mandatory): stimulation current
%       - location (optional): a logical matrix of the same size of atr and
%       ventr defining the region to stimulate. If location is not provided
%       the user is asked to select the regions to stimulate by mouse click.
%  In the case no both ExtStim and pkj fields are not passed to the
%  function, ExtStim will be automatically set up with a single stimulation
%  at the begining of the simulation with current equal to 50 and duration
%  equal to 1. The location will be selected by the user.
%
% - dt_visual: defines time interval for visualization. If dt_visual=Inf
% then no visualization is provided (default= 5 ms).
%
%- savefile: structure to specify if save the membrane potential to a
% binary file. The savefile structure must have the following subfields:
%       - filename: name of the file where the memebrane potential will be
%       saved
%       - dt: time interval for saving.
%
%- phase_xi_atr: xi parameter used for phase field computation in the 
% atrial geometry (see computePhase, default: 0.15 mm).
%
%- phase_xi_ventr: xi parameter used for phase field computation in the 
% ventricular geometry (see computePhase, default: 0.25 mm).
%
%- phase_dt: timestep used for phase field computation(default: 0.02 ms).
%
%- phase_thr: defines a phase field threshold below which voxels are not
%considered for computation of membrane potential (default: 1e-4).
%
%- save_state: structure to specify savings of the state of the
%simulation. The save_state structure may have two subfields:
%   - times: defines a vector of times at which saving the state
%     variables (default: [], no savings).
%   - filename: prefix name of mat files (default: 'sim').
%
%- load_state: filename of a mat file containing the state variables of the
%ionic models to use as initial condition.
%
%- ext_pot:structure to specify computation of extracellular potential by
%assuming infinite homogeneous conductor. The ext_pot structure may have
%the following subfiles:
%   - pos: a NeX3 matrix indicating the positions at which computing the
%   extracellular potential (mandatory).
%   - cond: diffusivity of the bath (default: 0.003 cm^2/ms)



if nargin<2
    options=[];
end

VoxelMat=atr | ventr;
nin_ventr=nnz(ventr);
nin_atr=nnz(atr);

if isfield(options,'domain')
    domain_ventr=domain(ventr);
    domain_atr=domain(atr);
else
    if isfield(options,'domain_ventr')
        domain_ventr=options.domain_ventr;
        if numel(domain_ventr)==numel(VoxelMat)
            domain_ventr=domain_ventr(ventr(:));
        end
    else
        domain_ventr=ones(nin_ventr,1);
    end

    if isfield(options,'domain_atr')
        domain_atr=options.domain_atr;
        if numel(domain_atr)==numel(VoxelMat)
            domain_atr=domain_atr(atr(:));
        end
    else
        domain_atr=ones(nin_atr,1);
    end
end


if isfield(options,'fib')
    fib_ventr=options.fib(ventr);
    fib_atr=options.fib(atr);
    plot_fib=options.fib; 
else
    if isfield(options,'fib_ventr')
        fib_ventr=options.fib_ventr;
        if numel(fib_ventr)==numel(VoxelMat)
            fib_ventr=fib_ventr(ventr);
        end
    else
        fib_ventr=false(nin_ventr,1);
    end
    
    if isfield(options,'fib_atr')
        fib_atr=options.fib_atr;
        if numel(fib_atr)==numel(VoxelMat)
            fib_atr=fib_atr(atr);
        end
    else
        fib_atr=false(nin_atr,1);
    end
    plot_fib=false(size(VoxelMat));
    plot_fib(atr)=fib_atr;
    plot_fib(ventr)=fib_ventr;
end


if isfield(options,'scalar')
    scalar_ventr=options.scalar(ventr);
    scalar_atr=options.scalar(atr);
else
    if isfield(options,'scalar_ventr')
        scalar_ventr=options.scalar_ventr;
        if numel(scalar_ventr)==numel(VoxelMat)
            scalar_ventr=scalar_ventr(ventr(:));
        end
    else
        scalar_ventr=zeros(nin_ventr,1);
    end

    if isfield(options,'scalar_atr')
        scalar_atr=options.scalar_atr;
        if numel(scalar_atr)==numel(VoxelMat)
            scalar_atr=scalar_atr(atr(:));
        end
    else
        scalar_atr=zeros(nin_atr,1);
    end
end

if isfield(options,'res')
    res=options.res;
else
    res=0.25;
end

if isfield(options,'fibers_ventr')
    f=options.fibers_ventr;
else
    f.x=ones(nin_ventr,1);
    f.y=zeros(nin_ventr,1);
    f.z=zeros(nin_ventr,1);
end

if isfield(options,'fibers_atr')
    fa=options.fibers_atr;
else
    fa.x=ones(nin_atr,1);
    fa.y=zeros(nin_atr,1);
    fa.z=zeros(nin_atr,1);
end

if isfield(options,'Dpar_ventr')
    Dpar_ventr=options.Dpar_ventr;
else
    Dpar_ventr=1.2e-3;
end

if isfield(options,'Dpar_atr')
    Dpar_atr=options.Dpar_atr;
else
    Dpar_atr=4.8e-3;
end

if isfield(options,'anisot_ratio_ventr')
    anisot_ratio_ventr=options.anisot_ratio_ventr;
else
    anisot_ratio_ventr=4;
end

if isfield(options,'anisot_ratio_atr')
    anisot_ratio_atr=options.anisot_ratio_atr;
else
    anisot_ratio_atr=8;
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
    if  isfield(options,'load_state')
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

if isfield(options,'phase_xi_ventr')
    xi_ventr=options.phase_xi_ventr/10;
else
    xi_ventr=0.025;
end

if isfield(options,'phase_xi_atr')
    xi_atr=options.phase_xi_atr/10;
else
    xi_atr=0.015;
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

if isfield(options,'model_ventr')
    model_ventr=options.model_ventr;
else
    model_ventr='Biasi_Tognetti';
end

if isfield(options,'model_atr')
    model_atr=options.model_atr;
else
    model_atr='Biasi_Tognetti';
end

%% Define simulation length
t=gpuArray(0:dt:T);
Nt=length(t);

%% Phase field and mask computation

phi_atr=computePhase(atr,res/10,phase_dt,xi_atr,true);
mask_phi_atr=phi_atr>thr_b;

phi_ventr=computePhase(ventr,res/10, phase_dt,xi_ventr,true);
mask_phi_ventr=phi_ventr>thr_b;

%% extrapolate fibers and domain and pixel intensity
[~,extInd_ventr]=computeSurface(ventr,res);
[~,extInd_atr]=computeSurface(atr,res);

fields=extrapField([f.x,f.y,f.z,domain_ventr,scalar_ventr,fib_ventr],ventr,extInd_ventr,mask_phi_ventr);
fx=fields(:,1);
fy=fields(:,2);
fz=fields(:,3);
fib_ventr=gpuArray(fields(:,6));
plot_fib2_ventr=false(size(VoxelMat));
plot_fib2_ventr(mask_phi_ventr)=logical(fib_ventr(mask_phi_ventr(:)));
ind_in_ventr= mask_phi_ventr & not(plot_fib2_ventr);
nin_ventr=nnz(ind_in_ventr);
domain_ventr=gpuArray(fields(ind_in_ventr,4));
scalar_ventr=gpuArray(fields(ind_in_ventr,5));


clear fields

fields_atr=extrapField([fa.x,fa.y,fa.z,domain_atr, scalar_atr, fib_atr],atr, extInd_atr, mask_phi_atr);
fax=fields_atr(:,1);
fay=fields_atr(:,2);
faz=fields_atr(:,3);
fib_atr=gpuArray(fields_atr(:,6));
plot_fib2_atr=false(size(VoxelMat));
plot_fib2_atr(mask_phi_atr)=logical(fib_atr(mask_phi_atr(:)));
ind_in_atr=mask_phi_atr & not(plot_fib2_atr);
nin_atr=nnz(ind_in_atr);
domain_atr=gpuArray(fields_atr(ind_in_atr,4));
scalar_atr=gpuArray(fields_atr(ind_in_atr,5));

clear fields_atr

%% elab purkinje system

if isfield(options,'pkj')%
    pkn_nodes=options.pkj.nodes;
    pkn_elem=options.pkj.elem;
    avn_x=pkn_nodes(1,1);
    avn_y=pkn_nodes(1,2);
    avn_z=pkn_nodes(1,3);
    avn_stim=false(size(VoxelMat));
    stimAvnSize=2;
    avn_stim(avn_y-stimAvnSize:avn_y+stimAvnSize,avn_x-stimAvnSize:avn_x+stimAvnSize,avn_z-stimAvnSize:avn_z+stimAvnSize)=true;
    avn_stim_map=avn_stim(ind_in_atr);
    term_mask=not(ismember(pkn_elem(:,2),pkn_elem(:,1)));
    term_ind=pkn_elem(term_mask,2);
    term_ind=[1;sort(term_ind)];
    distM=dist_pkj(pkn_nodes,pkn_elem);
    if isfield(options.pkj,'CV')
        pkj_vel=options.pkj.CV;
    else
        pkj_vel=340;
    end
    distM=distM*res/pkj_vel*100;
    if isfield(options.pkj,'rv_delay') || isfield(options.pkj,'lv_delay')
        if isfield(options.pkj,'rv_delay')
            rv_delay=options.pkj.rv_delay;
        else
            rv_delay=0;
        end
        if isfield(options.pkj,'lv_delay')
            lv_delay=options.pkj.lv_delay;
        else
            lv_delay=0;
        end
        distM(1,:)=distM(1,:)+[zeros(1,2) lv_delay*ones(1,(size(pkn_nodes,1)-2)/2) rv_delay*ones(1,(size(pkn_nodes,1)-2)/2)];
        distM(:,1)=distM(:,1)+[zeros(2,1); lv_delay*ones((size(pkn_nodes,1)-2)/2,1); rv_delay*ones((size(pkn_nodes,1)-2)/2,1)];
    end
    distM=distM(term_ind,term_ind);
    term_ind=term_ind(2:end);
    term_nodes=pkn_nodes(term_ind,:);
    term_map=round(sub2ind(size(VoxelMat),term_nodes(:,2),term_nodes(:,1),term_nodes(:,3)));
    tmp=zeros(size(VoxelMat));
    tmp(term_map)=2:length(term_nodes)+1;
    tmp=tmp(ind_in_ventr);
    if isfield(options.pkj,'fib_thr')
        fib_thr=options.pkj.fib_thr;
    else
        fib_thr=0.4;
    end
    tmp(scalar_ventr>fib_thr)=0;
    term_map2=find(tmp);
    good_term=tmp(term_map2);

    tmp=zeros(size(atr));
    avn_ind = sub2ind(size(atr),avn_y,avn_x,avn_z);
    tmp(avn_ind)=1;
    tmp=tmp(ind_in_atr);
    avn_map=find(tmp);
    plot_fib2_atr(avn_ind)=0; %remove fibrosis from avn point

    clear tmp

    if isfield(options.pkj,'stimTimes')
        pkn_stim=options.pkj.stimTimes;
    else
        pkn_stim=[];
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
    if isfield(options.pkj,'avn_anterograde_delay')
    avn_delayA=options.pkj.avn_anterograde_delay;
    else
        avn_delayA=80;
    end

    if isfield(options.pkj,'avn_retrograde_delay')
        avn_delayR=options.pkj.avn_retrograde_delay;
    else
        avn_delayR=80;
    end

    if isfield(options.pkj,'avn_stim_duration')
        avnStim_duration=options.pkj.avn_stim_duration;
    else
        avnStim_duration=2;
    end

    if isfield(options.pkj,'avn_current')
        avn_current=options.pkj.avn_current;
    else
        avn_current=200;
    end

    if isfield(options.pkj,'avn_erp')
        avn_erp=options.pkj.avn_erp;
    else
        avn_erp=300;
    end

    if isfield(options.pkj,'av_block')
        AVblock=options.pkj.av_block;
    else
        AVblock="none";
    end

else

    pkn_stim=[];
    warning('No Purkinje system provided, atria and ventricles are electrically disconnected!')
end
pkn_stim_n=round(pkn_stim/dt)+1;

%% system matrix
As_atr=gpuArray(heartMatrixFib(phi_atr,ind_in_atr,plot_fib2_atr,fax,fay,faz,Dpar_atr,anisot_ratio_atr,res));
As_ventr=gpuArray(heartMatrixFib(phi_ventr,ind_in_ventr,plot_fib2_ventr,fx,fy,fz,Dpar_ventr,anisot_ratio_ventr,res));
%% Set up external stimulation
for i=1:length(Stim)
    if Stim(i).toSelect
        disp(['Select stimulation region ', num2str(i), ' by mouse click']);
        Stim(i).ExtLoc=selectVolume(VoxelMat,res,5,false, plot_fib);
    end
    I_mask=zeros(size(VoxelMat));
    I_mask(Stim(i).ExtLoc & VoxelMat)=Stim(i).Iext;
    Stim(i).I_mask_atr=gpuArray(I_mask(ind_in_atr));
    Stim(i).I_mask_ventr=gpuArray(I_mask(ind_in_ventr));

    n_stim=length(Stim(i).ExtTimes);

    NT=round(Stim(i).ExtDuration/dt)+1;

    Stim(i).ext_stim=false(Nt,1);
    for j=1:n_stim
        ind=round(Stim(i).ExtTimes(j)/dt)+1;
        Stim(i).ext_stim(ind:ind+NT)=1;
    end
end
%% Initialization
[Vm_ventr0, state_model_ventr0]=feval(model_ventr);
nVar_ventr=length(state_model_ventr0);
[Vm_atr0, state_model_atr0]=feval(model_atr);
nVar_atr=length(state_model_atr0);
if isfield(options,'load_state')
    saved_state=load(options.load_state);
    Vm_ventr=saved_state.Vm_ventr;
    state_model_ventr=saved_state.state_model_ventr;
    if length(state_model_ventr)~=nVar_ventr
        error('number of loaded state variables do not correspond to the selected ventricular model');
    end
    Vm_atr=saved_state.Vm_atr;
    state_model_atr=saved_state.state_model_atr;
    if length(state_model_atr)~=nVar_atr
        error('number of loaded state variables do not correspond to the selected atrial model');
    end
else
    Vm_ventr=gpuArray(Vm_ventr0*ones(nin_ventr,1));
    state_model_ventr=cell(nVar_ventr,1);
    for ii=1:nVar_ventr
        state_model_ventr{ii}=state_model_ventr0{ii}*ones(nin_ventr,1);
    end
    Vm_atr=gpuArray(Vm_atr0*ones(nin_atr,1));
    state_model_atr=cell(nVar_atr,1);
    for ii=1:nVar_atr
        state_model_atr{ii}=state_model_atr0{ii}*ones(nin_atr,1);
    end
end

Ie_ventr=gpuArray(zeros(nin_ventr,1));
Ie_atr=gpuArray(zeros(nin_atr,1));

if ~isinf(visual)
    Vm_plot=NaN(size(VoxelMat));
    Vm_plot(ind_in_atr)=gather(Vm_atr);
    Vm_plot(ind_in_ventr)=gather(Vm_ventr);
    [FV1,extInd1]=computeSurface(plot_fib,res);
    hf=figure;
    hf.MenuBar='figure';
    [FV,extInd]=computeSurface(VoxelMat & ~plot_fib,res);
    h=PlotVoxel(FV,Vm_plot, extInd, parula,1,0);
    ttl=title('Time = 0 ms');
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
    fwrite(fid,gather(Vm_atr),'single');
    fwrite(fid,gather(Vm_ventr),'single');
    save(filename,'ind_in_atr','ind_in_ventr','plot_fib','atr','ventr','res','dt_save','t');
end
%% initializa PKJ
if  exist('term_map2','var')
    if isfield(options,'load_state')
        if isfield(saved_state,'state')
            state=saved_state.state;
            t_last=saved_state.t_last;
            Tact=saved_state.Tact;
            state_avn=saved_state.state_avn;
            t_last_avn=saved_state.t_last_avn;
            Tact_avn=saved_state.Tact_avn;
            if length(state)~=length(good_term)
                error(['Purkinje network in the state file does not '...
                    'correspond to Purkinje network specified in the ' ...
                    'options structure']);
            end
        else
            warning(['The loaded state file does not contain Purkinje states;'...
                ' Purkinje states will be initialized at rest'])
            state=gpuArray(zeros(size(good_term)));
            t_last=gpuArray(zeros(size(good_term)));
            Tact=Inf(size(good_term))';
            state_avn=0;
            t_last_avn=0;
            Tact_avn=Inf;
        end
    else
        state=gpuArray(zeros(size(good_term)));
        t_last=gpuArray(zeros(size(good_term)));
        Tact=Inf(size(good_term))';
        state_avn=0;
        t_last_avn=0;
        Tact_avn=Inf;
    end
else
    if isfield(options,'load_state')
        if isfield(saved_state,'state')
            warning(['Loaded state contains Purkinje network states but ' ...
                'no Purkinje network has been provided in the options' ...
                'structure. The simulation will be executed without ' ...
                'Purkinje network']);
        end
    end
end
    %% extracellular potential matrix
    if isfield(options,'ext_pot')
        pos = options.ext_pot.pos;
        if isfield(options.ext_pot,'cond')
            bath_cond_ventr=options.ext_pot.cond;
            bath_cond_atr=options.ext_pot.cond;
        else
            if isfield(options.ext_pot,'cond_ventr')
                bath_cond_ventr=options.ext_pot.cond_ventr;
            else
                bath_cond_ventr=0.0003;
            end

            if isfield(options.ext_pot,'cond_atr')
                bath_cond_atr=options.ext_pot.cond_atr;
            else
                bath_cond_atr=0.0012;
            end
        end
        LF_ventr=gpuArray(1/4/pi/bath_cond_ventr*LFvector(phi_ventr,mask_phi_ventr,plot_fib2_ventr,fx,fy,fz,Dpar_ventr,anisot_ratio_ventr,res,pos));
        LF_atr=gpuArray(1/4/pi/bath_cond_atr*LFvector(phi_atr,mask_phi_atr,plot_fib2_atr,fax,fay,faz,Dpar_atr,anisot_ratio_atr,res,pos));
        Vext=(NaN(Nt,size(pos,1)));
        Vext(1,:)=gather(LF_ventr*Vm_ventr+LF_atr*Vm_atr);
        figure;
        h_ext=plot(t,Vext);
        ax_ext=h_ext.Parent;
        ax_ext.XLim=[0 T];
    end


    %% RUN SIMULATION

    for i=1:Nt-1
        dv2dt_ventr=As_ventr*Vm_ventr;
        dv2dt_atr=As_atr*Vm_atr;
        for j=1:length(Stim)
            if Stim(j).ext_stim(i)
                dv2dt_ventr=dv2dt_ventr+Stim(j).I_mask_ventr;
                dv2dt_atr=dv2dt_atr+Stim(j).I_mask_atr;
            end
        end

        [Vm_ventr,state_model_ventr,dv2dt_ventr]=feval(model_ventr,Vm_ventr,state_model_ventr,Ie_ventr,dv2dt_ventr,domain_ventr,dt,scalar_ventr);
        [Vm_atr,state_model_atr,dv2dt_atr]=feval(model_atr,Vm_atr,state_model_atr,Ie_atr,dv2dt_atr,domain_atr,dt,scalar_atr);

        if toSave
            if mod(i,i_save)==0
                fwrite(fid,gather(Vm_atr),'single');
                fwrite(fid,gather(Vm_ventr),'single');
            end
        end

        if exist('Tact','var')
            active_avn=dv2dt_atr(avn_map)>0 & Vm_atr(avn_map)>-20;
            active=dv2dt_ventr(term_map2)>0 & Vm_ventr(term_map2)>-20;
            Tact_avn=min([Tact_avn; retro_delay+distM(good_term(active),1)],[],1);
            Tact_avn=Tact_avn-dt;        
            [sA,sV,state_avn,t_last_avn] = AVnodeFun(active_avn,Tact_avn<0,t(i),state_avn,t_last_avn,avn_delayA,avn_delayR,0,avnStim_duration,avn_erp);
            if Tact_avn<0
                Tact_avn=Inf;
            end

            if ismember(i,pkn_stim_n) || (sV && AVblock ~= "anterograde" && AVblock ~= "both")
                Tact=min(Tact,distM(1,good_term));
            end

            Tact=Tact-dt;
            [s_v,s_pkj,state,t_last]=PKJnodeFun(Tact'<0,active,t(i),state,t_last,ant_delay,retro_delay,pkj_duration,0,ERP);
            Tact(Tact<0)=Inf;
            Tact=min([Tact;distM(good_term(s_pkj),good_term)],[],1);
            Ie_ventr(:)=0;
            Ie_ventr(term_map2(s_v))=pkj_current;

            Ie_atr(:)=0;
            if(sA && AVblock ~= "retrograde" && AVblock ~= "both")
                Ie_atr(avn_stim_map)=avn_current;
            end
        end

        if ismember(dt*i,savestate_times)
            if exist('Tact','var')
                tmp=t_last;
                t_last=t_last-dt*i;
                tmp2=t_last_avn;
                t_last_avn=t_last_avn-dt*i;
                save([filename_state '_state_' num2str(dt*i)], 'Vm_ventr','state_model_ventr','Vm_atr','state_model_atr','Tact','t_last','state','Tact_avn','state_avn','t_last_avn');
                t_last=tmp;
                t_last_avn=tmp2;
            else
                save([filename_state '_state_' num2str(dt*i)], 'Vm_ventr','state_model_ventr','Vm_atr','state_model_atr');
            end
        end

        if exist('Vext','var')
            Vext(i+1,:)=LF_ventr*Vm_ventr+LF_atr*Vm_atr;
        end


        if mod(i,visual)==0
            Vm_plot(ind_in_ventr)=gather(Vm_ventr);
            Vm_plot(ind_in_atr)=gather(Vm_atr);
            h.CData=Vm_plot(extInd);
            ttl.String=['Time = ' num2str(dt*i) ' ms'];
            if exist('Vext','var')
                for pp=1:size(Vext,2)
                    h_ext(pp).YData=Vext(:,pp);
                end
            end
            drawnow
        end

        if exist('term_map2','var')
            check_pkj=nnz(~isinf(Tact))==0 && nnz(state<3 & state>0)==0 && nnz(pkn_stim(i:end))==0;
        else
            check_pkj=1;
        end

        if nnz(Vm_atr>=-40)<5 && nnz(Vm_ventr>=-40)<5 && check_pkj && state_avn~=1 && state_avn~=2
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

    if exist('Vext','var')
        if toSave
            save([options.savefile.filename '_ext'],'Vext','t');
        else
            save('ou_ext','Vext','t');
        end
    end

