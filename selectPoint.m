function [point_selected, voxel_ind]=selectPoint(VoxelMat,FV, extInd, res)
% selectPoint allows the user selection of point of the surface of the
% voxelized geometry by mouse click.
%
% point_selected=selectPoint(VoxelMat,FV, extInd, res) returns the position
% vector of the point selected.
% FV and extInd must be generated with the computeSurface function. The
% resolution res is mandatory.
%
% [point_selected, voxel_ind]=selectPoint(VoxelMat,FV, extInd, res) returns
% the coordinates of the surface voxel nearest to the selected point and
% the linear index associted with it.

[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
ind_surf=unique(extInd);
grx=gridx(ind_surf);
gry=gridy(ind_surf);
grz=gridz(ind_surf);

figure;
h=PlotVoxel(FV,VoxelMat,extInd,'b');
hold on
handle_scatter=scatter3([],[],[],'filled','o','MarkerFaceColor','red');
h.ButtonDownFcn = {@click,handle_scatter,1};

while true
    str=input('press:\n Q to exit \n L to add light \n','s');
    if str=='L'
        camlight;
        continue
    end
    if str=='Q'
        break;
    end
end
if nargout>1
    distances=(grx-handle_scatter.XData).^2+(gry-handle_scatter.YData).^2+(grz-handle_scatter.ZData).^2;
    [~,voxel_ind_in]=min(distances);
    voxel_ind=ind_surf(voxel_ind_in);
    point_selected(1)=grx(voxel_ind_in);
    point_selected(2)=gry(voxel_ind_in);
    point_selected(3)=grz(voxel_ind_in);
else
    point_selected=[handle_scatter.XData; handle_scatter.YData;handle_scatter.ZData];
end
end