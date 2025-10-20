function plotFibers(fx,fy,fz,VoxelMat, transmural,options)
% plotFibers(fx,fy,fz,VoxelMat, transmural,res) produces a line plot of
% epicardial and endocardial fiber orientation [fx,fy,fz]. Each line
% represents a local avergaed fiber orientation. The voxelized
% geometry VoxelMat and the transmural coordinate are needed.
% Additional options can be passd through the options structure:
% - res: resolution of the geometry (default=0.25 mm)
% - Nlines: number of lines to plot onto endo and epicardial surface
% (default=5000)
% - scale: length of lines (default=3)
% - lineWidth: width of lines (default=2)


if nargin<6
    options=[];
end

if isfield(options,'res')
    res=options.res;
else
    res=0.25;
end

if isfield(options,'Nlines')
    Nlines=options.Nlines;
else
    Nlines=5000;
end

if isfield(options,'scale')
    scale=options.scale;
else
    scale=3;
end

if isfield(options,'lineWidth')
    lineWidth=options.lineWidth;
else
    lineWidth=2;
end

if length(transmural)~=numel(VoxelMat) && length(transmural)==nnz(VoxelMat)
    tmp=nan(numel(VoxelMat),1);
    tmp(VoxelMat(:))=transmural;
    transmural =tmp;
end


[~,extInd]=computeSurface(VoxelMat,res);
[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
ind_surf=unique(extInd);
isSurf=false(numel(VoxelMat),1);
isSurf(ind_surf)=1;
fields=extrapField([fx, fy, fz],VoxelMat,extInd,10);
fx=reshape(fields(:,1),size(VoxelMat));
fy=reshape(fields(:,2),size(VoxelMat));
fz=reshape(fields(:,3),size(VoxelMat));
fx=imgaussfilt3(fx,5);
fy=imgaussfilt3(fy,5);
fz=imgaussfilt3(fz,5);
normf=reshape(vecnorm([fx(:),fy(:),fz(:)],2,2),size(VoxelMat));
fx_plot=fx./normf;
fy_plot=fy./normf;
fz_plot=fz./normf;
figure
hold on
[FV,extInd]=computeSurface(reshape(0.2<transmural & transmural<0.8 & VoxelMat(:),size(VoxelMat)),res);
PlotVoxel(FV,VoxelMat,extInd,[0.3 0.3 0.3],0.6,0);
ind_endo=transmural==0 & isSurf(:)==1;
ind_start=find(ind_endo);
ind_start=ind_start(randperm(length(ind_start),Nlines));
hq1=quiver3(gridx(ind_start),gridy(ind_start),gridz(ind_start),scale*fx_plot(ind_start),scale*fy_plot(ind_start),scale*fz_plot(ind_start),'r',...
    'LineWidth',lineWidth,'ShowArrowHead','off','MarkerEdgeColor','k');
hq1.AutoScale='off';

ind_epi=transmural==1 & isSurf(:)==1;
ind_start=find(ind_epi);
ind_start=ind_start(randperm(length(ind_start),Nlines));
hq2=quiver3(gridx(ind_start),gridy(ind_start),gridz(ind_start),scale*fx_plot(ind_start),scale*fy_plot(ind_start),scale*fz_plot(ind_start),'b',...
    'LineWidth',lineWidth,'ShowArrowHead','off','MarkerEdgeColor','k');
hq2.AutoScale='off';

end