function [FV_cut, extInd_cut,normal,d]=DefineCP(VoxelMat,res)
% [FV_cut, extInd_cut,normal,d]=DefineCP(VoxelMat,res) define a cutplane by
% selecting three points on the voxelized geometry. If more than 3 points
% are selected, only the last three are considered.
% FV_cut and extInd_cut can be used to plot scalar field on the cut voxelized
% geometries. normal and d define the normal to the cutplane and d is
% distance from the origin.

[FV,extInd]=computeSurface(VoxelMat,res);
disp('Select at least 3 points for cutplane definition')
points=selectMPoint(VoxelMat,FV,extInd,res);
p1=points(:,end-2);
p2=points(:,end-1);
p3=points(:,end);
normal = cross(p1 - p2, p1 - p3);
d = (p1(1)*normal(1) + p1(2)*normal(2) + p1(3)*normal(3));
[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
plane=normal(1)*gridx+normal(2)*gridy+normal(3)*gridz<d;
[FV_cut,extInd_cut]=computeSurface(VoxelMat & plane,res);
end