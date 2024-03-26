function [vol_handle, viewer]=VolSh(vol_col,VoxelMat,extInd,cmap,plot_fib)
% [vol_handle, viewer]=VolSh(vol_col,VoxelMat,extInd,cmap,plot_fib)
% VolSh produce a 3D volume plot of vol_col data on the geometry define by
% VoxelMat. vol_col can be a matrix of the same size of VoxelMat or a
% vector with dimension equal to non-zero entries of VoxelMat. extInd array as
% created by the computeSurface(VoxelMat,res) function is needed. Optional inputs are
% cmap (colormap, default=parula) and plot_fib(a matrix defining fibrotic
% voxels). If plot_fib is passed, black fibrosis is shown above the
% geometry (VoxelMat should be 0 where plot_fib is 1)

if nnz(size(vol_col,1,2,3)==size(VoxelMat))==3
    vol_col=vol_col(VoxelMat);
end
if nargin<4
    cmap=parula;
end
vol_col=reshape(extrapField(vol_col,VoxelMat,extInd,[],1),size(VoxelMat));
vol_col(isnan(vol_col))=0;
vol_handle=volshow(vol_col);
vol_handle.AlphaData=VoxelMat;
vol_handle.Colormap=cmap;
viewer=vol_handle.Parent;
viewer.BackgroundColor=[1 1 1];
viewer.BackgroundGradient='off';
if nargin>4
    vol_fib=volshow(plot_fib,Parent=viewer);
    vol_fib.Colormap=[0 0 0];
end