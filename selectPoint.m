function [point_selected, voxel_ind]=selectPoint(VoxelMat,FV, extInd, res)
% selectPoint allows the user selection of point of the
% voxelized geometry by mouse click.
%
% point_selected=selectPoint(VoxelMat,FV, extInd, res) returns the position
% vector of the point selected.
% FV and extInd must be generated with the computeSurface function. The
% resolution res is mandatory.
%
% [point_selected, voxel_ind]=selectPoint(VoxelMat,FV, extInd, res) returns
% the coordinates of the voxel nearest to the selected point and
% the linear index associted with it.

[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
ind_vol=find(VoxelMat);
grx=gridx(VoxelMat);
gry=gridy(VoxelMat);
grz=gridz(VoxelMat);

figure;
h=PlotVoxel(FV,VoxelMat,extInd,'b');
hold on
handle_scatter=scatter3([],[],[],'filled','o','MarkerFaceColor','red');
h.ButtonDownFcn = {@click,handle_scatter,1};
handle_scatter_int=scatter3([],[],[],'filled','o','MarkerFaceColor','green');

while true
    str=input(['press:\n Q to exit \n'...
        ' L to add light and update \n'...
        ' C to define a cutplane \n'...
        ' F to flip the current cutplane \n'...
        ' Z to remove the current cutplane \n'...
        'any other to update\n'],'s');
    if str=='L'
        camlight;
        continue
    end
    if str=='Q'
        break;
    end
    if str=='C'
        disp('Select at least 3 points for cutplane definition')
        h.ButtonDownFcn = {@click,handle_scatter_int,false};
        while true

            clc
            str_int=input('press:\n Q to exit \n L to add light and update \n B to delete last selected point \n any other to update\n','s');
            if str_int=='L'
                camlight;
                continue
            end
            if str_int=='Q'
                break;
            end

            if str_int=='B'
                if ~isempty(handle_scatter_int.XData)
                    handle_scatter_int.XData=handle_scatter_int.XData(1:end-1);
                    handle_scatter_int.YData=handle_scatter_int.YData(1:end-1);
                    handle_scatter_int.ZData=[handle_scatter_int.ZData(1:end-1)];
                end
                continue
            end
        end
        points=[handle_scatter_int.XData; handle_scatter_int.YData;handle_scatter_int.ZData];
        p1=points(:,end-2);
        p2=points(:,end-1);
        p3=points(:,end);
        normal = cross(p1 - p2, p1 - p3);
        d = (p1(1)*normal(1) + p1(2)*normal(2) + p1(3)*normal(3));
        plane=normal(1)*gridx+normal(2)*gridy+normal(3)*gridz<d;
        [FV_cut,~]=computeSurface(VoxelMat & plane,res);
        h.Faces=FV_cut.faces;
        h.Vertices=FV_cut.vertices;
        handle_scatter_int.XData=[];
        handle_scatter_int.YData=[];
        handle_scatter_int.ZData=[];
        h.ButtonDownFcn = {@click,handle_scatter,true};
        continue
    end

     if str=='F'
        if exist('plane','var')
            plane=not(plane);
            [FV_cut,~]=computeSurface(VoxelMat & plane,res);
            h.Faces=FV_cut.faces;
            h.Vertices=FV_cut.vertices;
        end
        continue
     end

     if str=='Z'
        if exist('plane','var')
            [FV_cut,~]=computeSurface(VoxelMat,res);
            clear plane
            h.Faces=FV_cut.faces;
            h.Vertices=FV_cut.vertices;
        end
        continue
    end

end
if nargout>1
    distances=(grx-handle_scatter.XData).^2+(gry-handle_scatter.YData).^2+(grz-handle_scatter.ZData).^2;
    [~,voxel_ind_in]=min(distances);
    voxel_ind=ind_vol(voxel_ind_in);
    point_selected(1)=grx(voxel_ind_in);
    point_selected(2)=gry(voxel_ind_in);
    point_selected(3)=grz(voxel_ind_in);
else
    point_selected=[handle_scatter.XData; handle_scatter.YData;handle_scatter.ZData];
end
end