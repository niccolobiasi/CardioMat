function [VoxelMat,vox_label]=vtk2vox(filename,res,labels)
% vtk2vox converts a vtk file (unstructured grid type) into a voxelized
% geometry.
% VoxelMat=vtk2vox(filename,res) returns the voxelized geometry VoxelMat
% generated from the vtk file (filename) with a resolution res.
% [VoxelMat, vox_label]=vtk2vox(filename,res) also returns the labels of
% the voxelized geometry as converted from the vtk file.
% [VoxelMat,vox_label]=vtk2vox(filename,res,labels) only converts the
% geometry domains identified with the labels reported in the array labels.
if nargin<2
    res=0.25;
end

if isempty(res)
    res=0.25;
end

VTK_DATA=read_vtk_3d(filename);
if nargin>2
    cells_in=ismember(VTK_DATA.tag, labels);
    vtk_elem=VTK_DATA.elements(cells_in,:);
    usedPoints=unique(vtk_elem(:));
    ind_points=NaN(length(VTK_DATA.points),1);
    ind_points(usedPoints)=1:length(usedPoints);
    vtk_points=VTK_DATA.points(usedPoints,:);
    vtk_elem=ind_points(vtk_elem);
else
    vtk_points=VTK_DATA.points;
    vtk_elem=VTK_DATA.elements;
end
bound=5;
Xmax=max(vtk_points(:,1))+bound;
Ymax=max(vtk_points(:,2))+bound;
Zmax=max(vtk_points(:,3))+bound;
Xmin=min(vtk_points(:,1))-bound;
Ymin=min(vtk_points(:,2))-bound;
Zmin=min(vtk_points(:,3))-bound;
x=Xmin:res:Xmax;
y=Ymin:res:Ymax;
z=Zmin:res:Zmax;
TR=triangulation(vtk_elem,vtk_points);
F=freeBoundary(TR);
VoxelMat=inpolyhedron(F,vtk_points,x,y,z,'tol',res/2);
seed=find(VoxelMat);
seed=seed(randi(length(seed)));
eik1=solveEikonal(VoxelMat,seed);
VoxelMat(eik1==-1 & VoxelMat(:))=0;
[FV,extInd]=computeSurface(VoxelMat,res);
if nargout>1
    if nargin>2
        vtk_label=VTK_DATA.tag(cells_in);
    else
        vtk_label=VTK_DATA.tag;
    end
    vox_label=nan(size(VoxelMat));
    [gridx,gridy,gridz]=meshgrid(x,y,z);
    queryPoints=[gridx(VoxelMat) gridy(VoxelMat) gridz(VoxelMat)];
    resElem=pointLocation(TR,queryPoints);
    resElem_plot=zeros(size(VoxelMat));
    resElem_plot(VoxelMat)=resElem;
    bad=isnan(resElem_plot);
    ind_bad=find(bad);
    [Ny, Nx, ~]=size(VoxelMat);
    while not(isempty(ind_bad))
        for ii=[1 -1 Ny -Ny Ny*Nx -Ny*Nx]
            upt_bad=not(bad(ind_bad+ii)) & VoxelMat(ind_bad+ii);
            upt=ind_bad(upt_bad);
            resElem_plot(upt)=resElem_plot(upt+ii);
            bad=isnan(resElem_plot);
            ind_bad=find(bad);
        end
    end
    vox_label(VoxelMat)=vtk_label(resElem_plot(VoxelMat));

    figure;
    PlotVoxel(FV,vox_label,extInd);
else
    figure;
    PlotVoxel(FV,VoxelMat,extInd,[0.8500 0.3250 0.0980]);
end
