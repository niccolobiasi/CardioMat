function [VoxelMat, int_scalar_mat]=adas2vox(adas_folder,adas_filename,res,shrink,bound,suppress_plot)
% adas2vox converts adas shells into voxelized geometry
%[VoxelMat, int_scalar_mat]=adas2vox(adas_folder,adas_filename)
%reads the ADAS .vtk shells with filename adas_filename in the folder
%adas_folder. The name of the vtk file should be in the format: 
% "adas_filename"_LV_ii_ENH_DE-MRI 2D.vtk
%where ii goes from 10 to 90.
%
%
%[VoxelMat, int_scalar_mat]=adas2vox(adas_folder,adas_filename,res) uses
%a specific resolution for the output of the voxelized geometry. The
%default value for the resolution is 0.25 mm.
%
%
%[VoxelMat, int_scalar_mat]=adas2vox(adas_folder,adas_filename,res, shrink) uses
%a specific shrink factor to be passed to the boundary MATLAB function. The
%default values for the shrink factor is 0.9.
%
%
%[VoxelMat,int_scalar_mat]=adas2vox(adas_folder,adas_filename,res,shrink,bound) uses a
%specific computational domain bound to delimit the cardiac domain. The
%default value for the bound is 5 mm.
%
%
%[VoxelMat,int_scalar_mat]=adas2vox(adas_folder,adas_filename,res,shrink,bound, true) 
% avoid plotting the ADAS normalized pixel intesity.
%
%
%This function employs the inpolyhedron function available in the MATLAB
%file exchange at:
% www.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume


%% Read VTK ADAS shell

Points=[];
Scalars=[];

for ii=1:9
    VTK_shell=read_vtk([adas_folder '/' adas_filename '_LV_' num2str(ii*10) '_ENH_DE-MRI 2D.vtk']);
    Points=[Points; VTK_shell.points];
    Scalars=[Scalars; VTK_shell.scalar];
end
%% Find a conforming 3D boundary of the points in the ADAS shells

if nargin<5
    bound=5;
end

if nargin<4
    shrink=0.9;
end

if nargin<3
    res=0.25;
end

bnd=boundary(Points,shrink);

Xmax=max(Points(:,1))+bound;
Ymax=max(Points(:,2))+bound;
Zmax=max(Points(:,3))+bound;

Xmin=min(Points(:,1))-bound;
Ymin=min(Points(:,2))-bound;
Zmin=min(Points(:,3))-bound;

x=Xmin:res:Xmax;
y=Ymin:res:Ymax;
z=Zmin:res:Zmax;

VoxelMat=inpolyhedron(bnd,Points,x,y,z);


%% interpolate scalar in ADAS shell into the voxelized geometry
if nargout>1
    F=scatteredInterpolant(Points,Scalars);
    [gridx,gridy,gridz]=meshgrid(x,y,z);
    int_scalar=F(gridx(VoxelMat),gridy(VoxelMat),gridz(VoxelMat));
    int_scalar_mat=NaN(size(gridx));
    int_scalar_mat(VoxelMat)=int_scalar;
end
%% remove unreacheable voxel
eik1=-1*ones(numel(VoxelMat),1);
k=1;
domain=zeros(size(VoxelMat));
bs=0;
while nnz(eik1(VoxelMat)<0)     
    seed=find(eik1==-1 & VoxelMat(:));
    if length(seed)<bs
        break;
    end
    seed=seed(randi(length(seed)));
    eik1=solveEikonal(VoxelMat,seed);
    sel=eik1>=0;
    Ns=nnz(sel);
    if Ns>bs
        bk=k;
        bs=Ns;
    end
    domain(sel)=k;
    k=k+1;
end
VoxelMat(domain~=bk)=0;
if k>1
    disp('The geometry is composed by multiple separated domains. Only the largest domain is mantained')
end

%% Plot normalized pixel intensity in the voxelized geometry

if nargin<6
    suppress_plot=false;
end

if not(suppress_plot)
    [FV, extInd]=computeSurface(VoxelMat,res,[Xmin,Ymin,Zmin]);
    figure;
    if nargout>1
        PlotVoxel(FV,int_scalar_mat,extInd,parula);
    else
        PlotVoxel(FV,VoxelMat,extInd,'r');
    end
    camlight
end
