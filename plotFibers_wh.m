function plotFibers_wh(fvx,fvy,fvz,fax,fay,faz,ventr,atr,transmural,phi,options)
% plotFibers_wh produces a line plot of whole-heart fiber orientation.
% Input arguments are ventricular fiber orientation (fvx,fvy,fvz), atrial
% fiber orinetation (fax,fay,faz), ventricular geometry (ventr), atrial
% geometry (atr), ventricular transmural coordinate (transmural), and
% atrial transmural coordinate (phi). An options structure can be passed
% for defining optional input arguments (see plotFiber help for details on
% the optional parameter).

if nargin<11
    options=[];
end

if ~isfield(options,'Nlines')
    options.Nlines=10000;
end

if length(transmural)==numel(ventr)
    transmural =transmural(ventr);
end

if length(phi)==numel(atr)
    phi =phi(atr);
end

VoxelMat=atr | ventr;
nin=nnz(VoxelMat);

fx=nan(nin,1);
fy=nan(nin,1);
fz=nan(nin,1);

fx(atr(VoxelMat))=fax;
fx(ventr(VoxelMat))=fvx;
fy(atr(VoxelMat))=fay;
fy(ventr(VoxelMat))=fvy;
fz(atr(VoxelMat))=faz;
fz(ventr(VoxelMat))=fvz;

phi_wh=nan(numel(VoxelMat),1);
phi_wh(atr(:))=phi;
phi_wh(ventr(:))=transmural;

plotFibers(fx,fy,fz,VoxelMat,phi_wh,options);
