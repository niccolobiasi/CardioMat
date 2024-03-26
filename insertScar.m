function plot_fib=insertScar(VoxelMat,scalar,res,suppress_plot)
% plot_fib=insertScar(VoxelMat,scalar) generates a logical matrix plot_fib
% ofthe same size of VoxelMat (i.e., the voxelized geometry). Each element
% of plot_fib equal to 1 corresponds to a scar/fibrotic element. The input
% parameter scalar determines the percentage of fibrotic/scar tissue. The
% scalar input can be a matrix of the same size of VoxelMat or a vector of
% size such that length(scalar)=nnz(VoxelMat). Each voxel in the voxelized
% geometry has a probability equal to the corresponding scalr value of
% being fibrotic/scar. Alternatively, the scalar input can be a single
% value, and in that case each voxel has the same probability of being
% fibrotic.
%
if nargin<4
    suppress_plot=false;
end

if numel(scalar)==numel(VoxelMat)
    scalar=scalar(VoxelMat);
elseif numel(scalar)==1
    scalar=scalar*ones(size(VoxelMat));
end

randNumber=rand(nnz(VoxelMat),1);
fib=randNumber<scalar;
plot_fib=false(size(VoxelMat));
plot_fib(VoxelMat)=fib;

if ~suppress_plot
    if nargin<3
        res=0.25;
    end
    [FV,extInd]=computeSurface(VoxelMat,res);
    figure;
    PlotVoxel(FV,double(plot_fib),extInd,gray,1,1);
end