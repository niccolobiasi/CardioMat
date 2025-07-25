function [point_selected, voxel_ind]=selectMPoint(VoxelMat,FV, extInd, res)
% selectMPoint allows the user selection of multiple points of the surface of the
% voxelized geometry by mouse click.
%
% point_selected=selectMPoint(VoxelMat,FV, extInd, res) returns the position
% vector of the points selected.
% FV and extInd must be generated with the computeSurface function. The
% resolution res is mandatory.
%
% [point_selected, voxel_ind]=selectMPoint(VoxelMat,FV, extInd, res) returns
% the coordinates of the surface voxel nearest to the selected points and
% the linear indices associted with it.

[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
ind_surf=unique(extInd);
grx=gridx(ind_surf);
gry=gridy(ind_surf);
grz=gridz(ind_surf);

hf=figure;
hf.MenuBar='figure';
h=PlotVoxel(FV,VoxelMat,extInd,'b');
hold on
handle_scatter=scatter3([],[],[],'filled','o','MarkerFaceColor','red');
replace_mode=false;
h.ButtonDownFcn = {@click,handle_scatter,replace_mode};

while true
    if replace_mode
        current_mode='single';
    else
        current_mode='multiple';
    end
    str=input(['press:\n Q to exit \n L to add light and update \n T to toggle point replacement (current mode: '  current_mode ')\n B to delete last selected point \n any other to update\n'],'s');
    if str=='L'
        camlight;
        continue
    end
    if str=='Q'
        break;
    end
    if str=='T'
        replace_mode=not(replace_mode);
        h.ButtonDownFcn = {@click,handle_scatter,replace_mode};
        continue
    end
    if str=='B'
        if ~isempty(handle_scatter.XData)
            handle_scatter.XData=handle_scatter.XData(1:end-1);
            handle_scatter.YData=handle_scatter.YData(1:end-1);
            handle_scatter.ZData=[handle_scatter.ZData(1:end-1)];
        end
        continue
    end
end
if nargout>1
    ntips=length(handle_scatter.XData);
    point_selected=zeros(3,ntips);
    for i=1:ntips
        distances=(grx-handle_scatter.XData(i)).^2+(gry-handle_scatter.YData(i)).^2+(grz-handle_scatter.ZData(i)).^2;
        [~,voxel_ind_in]=min(distances);
        voxel_ind=ind_surf(voxel_ind_in);
        point_selected(1,i)=grx(voxel_ind_in);
        point_selected(2,i)=gry(voxel_ind_in);
        point_selected(3,i)=grz(voxel_ind_in);
    end
else
    point_selected=[handle_scatter.XData; handle_scatter.YData;handle_scatter.ZData];
end
end