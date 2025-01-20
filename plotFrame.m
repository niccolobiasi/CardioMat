function h=plotFrame(filename_bin,filename_geom, cmap)
%plotFrame(filename_bin) read the binary file with name filename_bin and
%allows interactive visualzation of saved data. plotFrame search for a mat
%file in the same path and with the same name containing simulation
%metadata (geometry adn fibrosis). Otherwise the mat file filename can be
%passed as second argument to the plotFrame function.

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

load(filename_geom,"VoxelMat","ind_in",'plot_fib','res');
sizeVm=nnz(ind_in);
Vm_plot=nan(size(VoxelMat));
fseek(fid,0,'eof');
Nframe=ftell(fid)/4/sizeVm;
fseek(fid,0,'bof');


[FV,extInd]=computeSurface(VoxelMat & ~plot_fib,res);
FV_cut=FV;
extInd_cut=extInd;
[FV1,extInd1]=computeSurface(plot_fib,res);
FV_cut1=FV1;

% end


Vm=fread(fid,sizeVm,'single');
Vm_plot(ind_in)=Vm;


figure;
clf
h=PlotVoxel(FV_cut,Vm_plot, extInd_cut, cmap);
ax=h.Parent;
hold on
h1=PlotVoxel(FV_cut1,plot_fib,extInd1,'k');
handle_scatter_int=scatter3([],[],[],'filled','o','MarkerFaceColor','green');

% end
colorbar
clim([-85 50])
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
frame=1;

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
        ' Q to exit \n'],'s');
    if str=='Q'
        break;
    elseif str=='A'
        if frame==Nframe
            disp('Maximum frame achieved')
            continue;
        end
        Vm=fread(fid,sizeVm,'single');
        frame=frame+1;
    elseif str=='S'
        frame=input('Select frame:');
        if frame>Nframe || frame<1
            disp('Selcted frame does not exists')
            continue;
        end
        fseek(fid,4*sizeVm*(frame-1),'bof');
        Vm=fread(fid,sizeVm,'single');
    elseif str=='L'
        camlight;
        continue
    elseif str=='C'
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
        [Ny,Nx,Nz]=size(VoxelMat);
        x=(0.5:1:(Nx-0.5))*res;
        y=(0.5:1:(Ny-0.5))*res;
        z=(0.5:1:(Nz-0.5))*res;
        [gridx, gridy, gridz]=meshgrid(x,y,z);
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
    end


    Vm_plot(ind_in)=Vm;
    h.CData=Vm_plot(extInd_cut);
    title(['frame=' num2str(frame)])

end
fclose(fid);
end

