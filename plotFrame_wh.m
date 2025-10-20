function h=plotFrame_wh(filename_bin,filename_geom, cmap)
% plotFrame_wh(filename_bin) is the whole-heart equivalent of plotFrame.
% It reads the binary file with name filename_bin and allows interactive 
% visualzation of saved data. plotFrame_wh search for a mat file in the 
% same path and with the same name containing simulation metadata (geometry
% and fibrosis). Otherwise the mat file filename can be passed as second 
% argument to the plotFrame function.

fid=fopen(filename_bin,'r');

if nargin<2
    filename_geom=[filename_bin '.mat'];
end
if isempty(filename_geom)
    filename_geom=[filename_bin '.mat'];
end

if nargin<3
    cmap=parula;
end
if isempty(cmap)
    cmap=parula;
end

load(filename_geom,"ventr","atr","ind_in_atr",'ind_in_ventr','plot_fib','res','dt_save');
VoxelMat=ventr | atr;
sizeVm_atr=nnz(ind_in_atr);
sizeVm_ventr=nnz(ind_in_ventr);
Vm_plot=nan(size(VoxelMat));
fseek(fid,0,'eof');
Nframe=ftell(fid)/4/(sizeVm_atr+sizeVm_ventr);
fseek(fid,0,'bof');
t=0:dt_save:(Nframe-1)*dt_save;

[FV,extInd]=computeSurface(VoxelMat & ~plot_fib,res);
FV_cut=FV;
extInd_cut=extInd;
[FV1,extInd1]=computeSurface(plot_fib,res);
FV_cut1=FV1;


Vm=fread(fid,sizeVm_atr,'single');
Vm_plot(ind_in_atr)=Vm;
Vm=fread(fid,sizeVm_ventr,'single');
Vm_plot(ind_in_ventr)=Vm;


hf=figure;
hf.MenuBar='figure';
clf
h=PlotVoxel(FV_cut,Vm_plot, extInd_cut, cmap);
ax=h.Parent;
hold on
h1=PlotVoxel(FV_cut1,plot_fib,extInd1,'k',1,0);
handle_scatter_int=scatter3([],[],[],'filled','o','MarkerFaceColor','green');
handle_scatter_p=scatter3([],[],[],'filled','o','MarkerFaceColor','green');


[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
grx=[gridx(ind_in_atr); gridx(ind_in_ventr)];
gry=[gridy(ind_in_atr); gridy(ind_in_ventr)];
grz=[gridz(ind_in_atr); gridz(ind_in_ventr)];

colorbar
clim([-85 50])
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
frame=1;
ttl=title(['frame=' num2str(frame) '; time=' num2str(t(frame))]);

while true
    clc
    disp(['Number of saved frame: ' num2str(Nframe)])
    str=input(['press:\n A to advance a frame \n' ...
        ' S to select a frame \n' ...
        ' L to add light \n'...
        ' C to define a cutplane \n'...
        ' F to flip the current cutplane \n'...
        ' Z to remove the current cutplane \n'...
        ' T to toogle fibrosis visibility \n'...
        ' M to change colormap limits \n'...
        ' P to plot action potential \n'...
        ' Q to exit \n'],'s');
    if str=='Q'
        break;
    elseif str=='A'
        if frame==Nframe
            disp('Maximum frame achieved')
            continue;
        end
        Vm_atr=fread(fid,sizeVm_atr,'single');
        Vm_ventr=fread(fid,sizeVm_ventr,'single');
        frame=frame+1;
    elseif str=='S'
        frame=input('Select frame:');
        if frame>Nframe || frame<1
            disp('Selcted frame does not exists')
            continue;
        end
        fseek(fid,4*(sizeVm_atr+sizeVm_ventr)*(frame-1),'bof');
        Vm_atr=fread(fid,sizeVm_atr,'single');
        Vm_ventr=fread(fid,sizeVm_ventr,'single');
    elseif str=='L'
        camlight;
        continue
    elseif str=='C'
        h.ButtonDownFcn = {@click,handle_scatter_int,false};
        while true

            clc
            disp('Select at least 3 points for cutplane definition')
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
        if h1.Visible
            [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane,res);
        else
            [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane & ~plot_fib,res);
        end
        [FV_cut1,~]=computeSurface(plot_fib & plane,res);
        h.Faces=FV_cut.faces;
        h.Vertices=FV_cut.vertices;
        h1.Faces=FV_cut1.faces;
        h1.Vertices=FV_cut1.vertices;
        h.CData=Vm_plot(extInd_cut);
        handle_scatter_int.XData=[];
        handle_scatter_int.YData=[];
        handle_scatter_int.ZData=[];
        h.ButtonDownFcn =[];
        continue
    elseif str=='F'
        if exist('plane','var')
            plane=not(plane);
            if h1.Visible
                [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane,res);
            else
                [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane & ~plot_fib,res);
            end
            [FV_cut1,~]=computeSurface(plot_fib & plane,res);
            h.Faces=FV_cut.faces;
            h.Vertices=FV_cut.vertices;
            h1.Faces=FV_cut1.faces;
            h1.Vertices=FV_cut1.vertices;
            h.CData=Vm_plot(extInd_cut);
        end
        continue
    elseif str=='Z'
        if exist('plane','var')
            if h1.Visible
                FV_cut=FV;
                extInd_cut=extInd;
            else
                [FV_cut,extInd_cut]=computeSurface(VoxelMat & ~plot_fib,res);
            end
            FV_cut1=FV1;
            clear plane
            h.Faces=FV_cut.faces;
            h.Vertices=FV_cut.vertices;
            h1.Faces=FV_cut1.faces;
            h1.Vertices=FV_cut1.vertices;
            h.CData=Vm_plot(extInd_cut);
        end
        continue
    elseif str=='T'
        h1.Visible=not(h1.Visible);
        if h1.Visible
            if exist('plane','var')
                [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane,res);
                [FV_cut1,~]=computeSurface(plot_fib & plane,res);
            else
                FV_cut=FV;
                extInd_cut=extInd;
                FV_cut1=FV1;
            end
        else
            if exist('plane','var')
                [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane & ~plot_fib,res);
            else
                [FV_cut,extInd_cut]=computeSurface(VoxelMat & ~plot_fib,res);
            end
        end
        h.Faces=FV_cut.faces;
        h.Vertices=FV_cut.vertices;
        h.CData=Vm_plot(extInd_cut);
        h1.Faces=FV_cut1.faces;
        h1.Vertices=FV_cut1.vertices;
        continue
    elseif str=='M'
        clc
        strLim=input('Write colormap limits between square brackets \n');
        ax.CLim=strLim;
        continue





    elseif str=='P'
        h.ButtonDownFcn = {@click,handle_scatter_p,true};
        figure;
        ax=gca;
        ax.XLim=[0 t(end)];
        ax.YLim=[-90 50];
        while true
            clc
            disp('Select point for action potential plot')
            disp(['Current mode: ' ax.NextPlot])
            str_int=input(['press:\n Q to exit ...' ...
                '\n H to hold previous plot' ...
                '\n R to replace previous plots'...
                '\n any other to update\n'],'s');
            if str_int=='Q'
                break;
            elseif str_int=='H'
                ax.NextPlot='add';
                continue
            elseif str_int=='R'
                ax.NextPlot='replacechildren';
                continue
            end

            distances=(grx-handle_scatter_p.XData(end)).^2+(gry-handle_scatter_p.YData(end)).^2+(grz-handle_scatter_p.ZData(end)).^2;
            [~,voxel_ind_in]=min(distances);
            Vm_array=NaN(Nframe,1);
            for k=1:Nframe
                fseek(fid,4*(sizeVm_atr+sizeVm_ventr)*(k-1)+4*(voxel_ind_in-1),'bof');
                Vm_array(k)=fread(fid,1,'single');
            end

            plot(ax,t,Vm_array,'LineWidth',1.5);

        end
        handle_scatter_p.XData=[];
        handle_scatter_p.YData=[];
        handle_scatter_p.ZData=[];
        h.ButtonDownFcn =[];
        continue

    end

    Vm_plot(ind_in_atr)=Vm_atr;
    Vm_plot(ind_in_ventr)=Vm_ventr;
    h.CData=Vm_plot(extInd_cut);
    ttl.String=['frame=' num2str(frame) '; time=' num2str(t(frame))];

end
fclose(fid);
end

