function [nodes,elem,act_times]=createPurkinje_biv(VoxelMat,transmural,apicobasal,intraventricular, options)
% [nodes,elem,act_times]=createPurkinje_biv(VoxelMat,transmural,apicobasal,intraventricular, options)
% generates a physiological biventricular Purkinje network for the voxelized geometry VoxelMat. 
%
% This function calls the createPurkinje function two times (one for each ventricle).
% The intraventricular coordinate is needed to separate the two ventricles.
% The two generated tree are then connected to a user-selected common node.
% The possible options are the same of the createPurkinje function and are
% passed to the createPurkinje 

if nargin<5
    options=[];
end

if isfield(options,'plot')
    doplot=options.plot;
else
    doplot=true;
end
if isfield(options,'res')
    res=options.res;
else
    res=0.25;
end
if isfield(options,'cond_vel')
    cond_vel=options.cond_vel;
else
    cond_vel=340;
end

options.plot=false;
intraventricular=reshape(intraventricular,size(VoxelMat));
LV=VoxelMat & intraventricular<=0.5;
disp('Generating left ventricle Purkinje network...')
[lv_nodes,lv_elem,lv_act_times]=createPurkinje(LV,transmural,apicobasal,options);

RV=VoxelMat & intraventricular>0.5;
disp('Generating right ventricle Purkinje network...')
[rv_nodes,rv_elem,rv_act_times]=createPurkinje(RV,transmural,apicobasal,options);

[FV, extInd]=computeSurface(VoxelMat,1);
disp('Select His starting point')
[root, ~]=selectPoint(VoxelMat,FV, extInd, 1);
nodes=[root;lv_nodes; rv_nodes];
root_lv=find(lv_act_times==0)+1;
root_rv=find(rv_act_times==0)+1+length(lv_nodes);
elem=[1 root_lv; 1 root_rv; lv_elem+1;rv_elem+length(lv_elem)+2];
pathLen_new=vecnorm(nodes([1 1],:)-nodes([root_lv root_rv],:),2,2);
act_times_new=pathLen_new*res/cond_vel*100;
act_times=[0;lv_act_times+act_times_new(1);rv_act_times+act_times_new(2)];
if doplot
    plotPurkinje(VoxelMat,transmural,nodes,elem,act_times);
end