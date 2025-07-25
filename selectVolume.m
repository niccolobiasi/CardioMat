function plot_selected=selectVolume(VoxelMat,res,R,replace_mode, plot_fib)
% selectVolume allows the user selection of portion of volume in a
% voxelized geometry by mouse click.
% plot_selected=selectSurface(VoxelMat,res,R) returns a matrix
% plot_selected of the same size of VoxelMat (the voxelized geometry). FV
% and extInd are generated inside the function by calling the
% computeSurface function. The resolution res is mandatory. The parameter R
% is mandatory and is used as the radius of the selection sphere around the
% user-selected points. The parameter replace_mode is optional and if set 
% true initialize the selection as single point (replacing the old one).
% Optionally, the user can specify a plot_fib matrix to define
% fibrotic/scar voxels in the geometry that should not be considered for the
% selection
if nargin<4
    replace_mode=false;
end

if nargin<5
    plot_fib=false(size(VoxelMat));
end


[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
ind_sel=VoxelMat & not(plot_fib);
grx=gridx(ind_sel);
gry=gridy(ind_sel);
grz=gridz(ind_sel);


R2=R^2;

selected_points =false(nnz(ind_sel),1);
plot_selected=false(size(VoxelMat));

[FV_sel,extInd_sel]=computeSurface(ind_sel,res);
[FV_fib,extInd_fib]=computeSurface(plot_fib & VoxelMat,res);
hf=figure;
hf.MenuBar='figure';
h=PlotVoxel(FV_sel,plot_selected,extInd_sel);
extInd_sel_cut=extInd_sel;
clim([0 1])
hold on
handle_scatter=scatter3([],[],[],'filled','o','MarkerFaceColor','red');
handle_scatter_int=scatter3([],[],[],'filled','o','MarkerFaceColor','green');
% if ~isempty(extInd_fib)
    h1=PlotVoxel(FV_fib,plot_fib,extInd_fib,'k',1,0);
% end
h.ButtonDownFcn = {@click,handle_scatter,replace_mode};

while true
    prev=selected_points;
    if replace_mode
        current_mode='single';
    else
        current_mode='multiple';
    end
    str=input(['press:\n Q to exit \n'...
        ' L to add light and update \n'...
        ' R to change the selection sphere size (current R value: ' num2str(R) ')\n'...
        ' T to toggle point replacement (current mode: '  current_mode ')\n' ...
        ' B to delete last selected point \n'...
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
        [FV_sel_cut,extInd_sel_cut]=computeSurface(ind_sel & plane,res);
        h.Faces=FV_sel_cut.faces;
        h.Vertices=FV_sel_cut.vertices;
        h.CData=double(plot_selected(extInd_sel_cut));
        % if ~isempty(extInd_fib)
            [FV_fib_cut,~]=computeSurface(plot_fib & plane,res);
            h1.Faces=FV_fib_cut.faces;
            h1.Vertices=FV_fib_cut.vertices;
        % end
        handle_scatter_int.XData=[];
        handle_scatter_int.YData=[];
        handle_scatter_int.ZData=[];
        h.ButtonDownFcn = {@click,handle_scatter,replace_mode};
        continue
    end

    if str=='F'
        if exist('plane','var')
            plane=not(plane);
            [FV_sel_cut,extInd_sel_cut]=computeSurface(ind_sel & plane,res);
            h.Faces=FV_sel_cut.faces;
            h.Vertices=FV_sel_cut.vertices;
            h.CData=double(plot_selected(extInd_sel_cut));
            [FV_fib_cut,~]=computeSurface(plot_fib & plane,res);
            h1.Faces=FV_fib_cut.faces;
            h1.Vertices=FV_fib_cut.vertices;
        end
        continue
    end

    if str=='Z'
        if exist('plane','var')
            [FV_sel_cut,extInd_sel_cut]=computeSurface(ind_sel,res);
            clear plane
            h.Faces=FV_sel_cut.faces;
            h.Vertices=FV_sel_cut.vertices;
            h.CData=double(plot_selected(extInd_sel_cut));
            [FV_fib_cut,~]=computeSurface(plot_fib,res);
            h1.Faces=FV_fib_cut.faces;
            h1.Vertices=FV_fib_cut.vertices;
        end
        continue
    end


    ntips=length(handle_scatter.XData);
    for i=1:ntips
        new_sel=(grx-handle_scatter.XData(i)).^2+(gry-handle_scatter.YData(i)).^2+(grz-handle_scatter.ZData(i)).^2<R2;
        selected_points=or(selected_points, new_sel);
    end
    plot_selected(ind_sel)=selected_points;
    handle_scatter.XData=[];
    handle_scatter.YData=[];
    handle_scatter.ZData=[];
    h.CData=double(plot_selected(extInd_sel_cut));
    if ntips>0
        str2=input('Confirm volume (D to delete)?','s');
        if str2=='D'
            selected_points=prev;
            plot_selected(ind_sel)=selected_points;
            h.CData=double(plot_selected(extInd_sel_cut));
        end
    end




end