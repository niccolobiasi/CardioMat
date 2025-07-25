function PlotVoxelE(VoxelMat,Vm_plot,res, plot_fib,cmap,alpha,addLight)
% PlotVoxelE(VoxelMat, Vm_plot, res) allows interactive visualization of
% the voxelized 3D image Vm_plot defined only on the geometry VoxelMat.
%
% VoxelMat is the 3D voxelized geometry.
% Vm_plot can be a matrix of the same size of VoxelMat or a vector with
% length equal to nnz(VoxelMat).
% res is the resolution.
%
% PlotVoxelE(VoxelMat,Vm_plot,res, plot_fib,cmap,alpha,addLight) allows to
% additionally define a plot_fib matrix for fibrosis representation, the
% colormap cmap (default: parula), the alpha value for surface voxels 
% (default:1), and if light should be added or not (default:true)
%
%

hf=figure;
hf.MenuBar='figure';

if size(Vm_plot,2)==1
    tmp=Vm_plot;
    Vm_plot=nan(size(VoxelMat));
    Vm_plot(VoxelMat)=tmp;
    clear tmp
end
    

if nargin<4
    plot_fib=false(size(VoxelMat));
end

if nargin<7
    addLight=true;
end



if nargin<5
    cmap=parula(256);
end

VoxelMat=VoxelMat & ~plot_fib;

[FV,extInd]=computeSurface(VoxelMat,res);
[FV1,extInd1]=computeSurface(plot_fib,res);


if nargin<6
    alpha=1;
elseif size(alpha,1)==size(Vm_plot,1) && size(alpha,2)==size(Vm_plot,2)
    alpha=alpha(extInd)';
elseif size(alpha,1)==numel(Vm_plot)
    alpha=alpha(extInd);
end

h=PlotVoxel(FV,Vm_plot,extInd,cmap,alpha,addLight);
ax=h.Parent;
hold on
h1=PlotVoxel(FV1,plot_fib,extInd1,'k',1,0);





handle_scatter_int=scatter3([],[],[],'filled','o','MarkerFaceColor','green');
while true
    clc
    str=input(['press: \n'...
        ' L to add light \n'...
        ' C to define a cutplane \n'...
        ' F to flip the current cutplane \n'...
        ' Z to remove the current cutplane \n'...
        ' T to toogle fibrosis visibility \n'...
                ' M to change colormap limits \n'...
        ' Q to exit \n'],'s');
    if str=='Q'
        break;
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
        [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane,res);
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
            [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane,res);
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
            FV_cut=FV;
            extInd_cut=extInd;
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
        continue
    elseif str=='M'
        strLim=input('write max min in square brackets');
        ax.CLim=strLim;
    end


end