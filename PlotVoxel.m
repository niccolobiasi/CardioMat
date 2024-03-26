function vol_handle=PlotVoxel(FV,Vmplot,extInd,cmap,alpha,addLight)
% PlotVoxel((FV,Vmplot,myExtInd) generates a surface plot of a voxelized 
% geometry. The FV structure and extInd array must be generated with the
% computeSurface function.
%
% vol_handle=PlotVoxel(FV,Vmplot,extInd) also returns the handle to the
% surface plot.
%
% vol_handle=PlotVoxel(FV,Vmplot,extInd,cmap) allows to specify a colormap
% cmap. the default colormap is parula(256). cmap can also be a char
% indicating a color to be used for all the surfaces with the default
% matlab conventions (e.g., 'r' is red)
%
% vol_handle=PlotVoxel(FV,Vmplot,extInd,cmap,alpha) allows to specify an
% alpha value (trasparency) for each voxel. If a scalar value is provided
% the same value is used for all the voxels. The default value for alpha is
% 1. Note that the alpha value only applies to surface voxels. Internal
% voxels are not rendered (use VolSh at this scope).
%
%
%vol_handle=PlotVoxel(FV,Vmplot,extInd,cmap,alpha,addLight) just allows to
%add light in the plot by calling camlight matlab function
%
%



if nargin<6
    addLight=true;
end

if nargin<5
    alpha=1;
elseif size(alpha,1)==size(Vmplot,1) && size(alpha,2)==size(Vmplot,2)
    alpha=alpha(extInd)';
elseif size(alpha,1)==numel(Vmplot)
    alpha=alpha(extInd);
end

if nargin<4
    cmap=parula(256);
end

if islogical(Vmplot)
    Vmplot=double(Vmplot);
end

vol_handle=patch("Faces",FV.faces,'Vertices',FV.vertices,...
    'CData',Vmplot(extInd),'EdgeColor','none','FaceColor','flat', ...
    'DiffuseStrength', 0.6,'backfacelighting',....
    'reverselit','ambientstrength', 0.4,...
    'FaceLighting', 'flat',...
    'SpecularStrength', 0.4,...
    'SpecularExponent', 10,...
    'SpecularColorReflectance', 1.0, ...
    'EdgeLighting','none');
if ischar(cmap)
    vol_handle.FaceColor=cmap;
elseif size(cmap,1)==1
    vol_handle.FaceColor=cmap;
else
    colormap(cmap)
    colorbar
end
view(3);
daspect([1,1,1]);

vol_handle.FaceVertexAlphaData=alpha;

vol_handle.FaceAlpha='flat';
vol_handle.AlphaDataMapping='none';

if addLight==1
    camlight
end
cameratoolbar('SetCoordSys','none');
cameratoolbar('show');
end