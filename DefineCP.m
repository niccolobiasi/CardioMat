function [plane, FV_cut, extInd_cut,normal,d]=DefineCP(VoxelMat,res,plot_fib)
% [FV_cut, extInd_cut,normal,d]=DefineCP(VoxelMat,res) define a cutplane by
% selecting three points on the voxelized geometry. If more than 3 points
% are selected, only the last three are considered.
% FV_cut and extInd_cut can be used to plot scalar field on the cut voxelized
% geometries. normal and d define the normal to the cutplane and d is
% distance from the origin.
if nargin<3
    plot_fib=false(size(VoxelMat));
end

[FV,extInd]=computeSurface(VoxelMat,res);
[FV1,extInd1]=computeSurface(plot_fib,res);

[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);

figure;
h=PlotVoxel(FV,VoxelMat,extInd,'b');
hold on
h1=PlotVoxel(FV1,plot_fib,extInd1,'k');
handle_scatter=scatter3([],[],[],'filled','o','MarkerFaceColor','red');
replace_mode=false;
h.ButtonDownFcn = {@click,handle_scatter,replace_mode};

while true
    if replace_mode
        current_mode='single';
    else
        current_mode='multiple';
    end
    clc
    disp('Select at least 3 points for cutplane definition')
    str=input(['press:\n L to add light and update \n'...
        ' T to toggle point replacement (current mode: '  current_mode ')\n'...
        ' B to delete last selected point \n'...
        ' any other to update \n '],'s');
    if str=='L'
        camlight;
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
    points=[handle_scatter.XData; handle_scatter.YData;handle_scatter.ZData];
    p1=points(:,end-2);
    p2=points(:,end-1);
    p3=points(:,end);
    normal = cross(p1 - p2, p1 - p3);
    d = (p1(1)*normal(1) + p1(2)*normal(2) + p1(3)*normal(3));
    plane=normal(1)*gridx+normal(2)*gridy+normal(3)*gridz<d;
    [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane,res);
    [FV1_cut,~]=computeSurface(plot_fib & plane,res);
    h.Faces=FV_cut.faces;
    h.Vertices=FV_cut.vertices;
    h1.Faces=FV1_cut.faces;
    h1.Vertices=FV1_cut.vertices;
    clc
    str2=input('Press F to flip the cutplane any other to continue \n','s');
    if str2=='F'
        plane=not(plane);
        [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane,res);
        [FV_cut1,~]=computeSurface(plot_fib & plane,res);
        h.Faces=FV_cut.faces;
        h.Vertices=FV_cut.vertices;
        h1.Faces=FV_cut1.faces;
        h1.Vertices=FV_cut1.vertices;
    end
    clc
    str3=input([' Press D to delete cutplane and repeat points selection' ...
        '\n press any other to confirm \n '],'s');
    if str3=='D'
        h.Faces=FV.faces;
        h.Vertices=FV.vertices;
        h1.Faces=FV1.faces;
        h1.Vertices=FV1.vertices;
        continue;
    else
        break;
    end

end