function rotational=computeRotational(VoxelMat,apicobasal,doplot)
% rotational=computeRotational(VoxelMat,apicobasal,doplot) computes the
% rotational coordinate (nnz(VoxelMat)X1) in the ventricle described by
% VoxelMat. The apicobasal coordinate is required as a way to define apex
% and base. 
% rotational=computeRotational(VoxelMat,apicobasal,true) also plot the
% rotational coordinate.
% Note that the user is asked to select a septal point in which the
% rotational coordinate will be 0.
if nargin<3
    doplot=false;
end

[FV,extInd]=computeSurface(VoxelMat,1);


[ii,jj,kk]=ind2sub(size(VoxelMat),find(apicobasal==0));
apex=mean([jj ii kk]);
[ii,jj,kk]=ind2sub(size(VoxelMat),find(apicobasal==1));
base=mean([jj ii kk]);
disp('Select a septal point (rotational=0)')
septal=selectPoint(VoxelMat,FV,extInd,1)';
normal = cross(base - apex, base - septal);
d = (base(1)*normal(1) + base(2)*normal(2) + base(3)*normal(3));

[ii, jj, kk]=ind2sub(size(VoxelMat),find(VoxelMat));
front=(normal(1)*jj+normal(2)*ii+normal(3)*kk-d)>=0;
plot_front=false(size(VoxelMat));
plot_front(VoxelMat)=front;
plot_back=~plot_front & VoxelMat;

normal2=cross(normal,base-apex);
d2=(base(1)*normal2(1) + base(2)*normal2(2) + base(3)*normal2(3));

[~,extInd_front]=computeSurface(plot_front,1);
[~,extInd_back]=computeSurface(plot_back,1);
ind_surf=unique(extInd);
ind_front=unique(extInd_front);
ind_back=unique(extInd_back);
isSurf=false(size(VoxelMat));
isSurf_front=false(size(VoxelMat));
isSurf_back=false(size(VoxelMat));
isSurf(ind_surf)=1;
isSurf_front(ind_front)=1;
isSurf_back(ind_back)=1;

newSfront=isSurf_front & ~isSurf;
ind_new_front=find(newSfront);
septalS=(normal2(1)*jj(newSfront(VoxelMat))+normal2(2)*ii(newSfront(VoxelMat))+normal2(3)*kk(newSfront(VoxelMat))-d2)>=0;
seed=ind_new_front(septalS);
eik1=solveEikonal(plot_front,seed);
seed=ind_new_front(~septalS);
eik2=solveEikonal(plot_front,seed);
rot1=pi*eik2./(eik1+eik2);


newSback=isSurf_back & ~isSurf;
ind_new_back=find(newSback);
septalS=(normal2(1)*jj(newSback(VoxelMat))+normal2(2)*ii(newSback(VoxelMat))+normal2(3)*kk(newSback(VoxelMat))-d2)<0;
seed=ind_new_back(septalS);
eik1=solveEikonal(plot_back,seed);
seed=ind_new_back(~septalS);
eik2=solveEikonal(plot_back,seed);
rot2=-pi*eik1./(eik1+eik2);

rotational=zeros(nnz(VoxelMat),1);
rotational(front)=rot1(plot_front);
rotational(~front)=rot2(plot_back);
if doplot
    figure;
    rota=nan(size(VoxelMat));
    rota(VoxelMat)=rotational;
    PlotVoxel(FV,rota,extInd,hsv);
end
