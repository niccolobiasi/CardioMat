function plot_selected=selectSurface(VoxelMat,FV, extInd, res,R,replace_mode)
% selectSurface allows the user selection of portion of surfaces in a
% voxelized geometry by mouse click.
% plot_selected=selectSurface(VoxelMat,FV, extInd, res,R) returns a matrix
% plot_selected of the same size of VoxelMat (the voxelized geometry). FV
% and extInd must be generated with the computeSurface function. The
% resolution res is mandatory. The parameter R is mandatory and is used as
% the radius of the selection sphere around the user-selected points.
% The parameter replace_mode is optional and if set true initialize the
% selection as single point (replacing the old one). 

if nargin<6
    replace_mode=false;
end

[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
ind_surf=unique(extInd);
grx=gridx(ind_surf);
gry=gridy(ind_surf);
grz=gridz(ind_surf);

R2=R^2;

selected_points =false(length(ind_surf),1);
plot_selected=false(size(VoxelMat));
figure;
h=PlotVoxel(FV,plot_selected,extInd);
clim([0 1])
hold on
handle_scatter=scatter3([],[],[],'filled','o','MarkerFaceColor','red');
h.ButtonDownFcn = {@click,handle_scatter,replace_mode};

while true
    prev=selected_points;
    if replace_mode
        current_mode='single';
    else
        current_mode='multiple';
    end
    str=input(['press:\n Q to exit \n L to add light and update \n R to change the selection sphere size (current R value: ' num2str(R) ')\n T to toggle point replacement (current mode: '  current_mode ')\n B to delete last selected point \n any other to update\n'],'s');
    if str=='L'
        camlight;
        continue
    end
    if str=='Q'
        break;
    end
    if str=='R'
        R=input('new R value:');
        R2=R^2;
        continue
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


    ntips=length(handle_scatter.XData);
    for i=1:ntips
        new_sel=(grx-handle_scatter.XData(i)).^2+(gry-handle_scatter.YData(i)).^2+(grz-handle_scatter.ZData(i)).^2<R2;
        selected_points=or(selected_points, new_sel);
    end
    plot_selected(ind_surf)=selected_points;
    handle_scatter.XData=[];
    handle_scatter.YData=[];
    handle_scatter.ZData=[];
    h.CData=double(plot_selected(extInd));
    if ntips>0
        str2=input('Confirm surface (D to delete)?','s');
        if str2=='D'
            selected_points=prev;
            plot_selected(ind_surf)=selected_points;
            h.CData=double(plot_selected(extInd));
        end
    end




end