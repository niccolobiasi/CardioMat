function ablation=VirtualAblation(VoxelMat,int_scalar_mat,res,transmural,adas_folder,adas_filename)
% ablation=VirtualAblation(VoxelMat,int_scalar_mat,res,transmural) allows
% the user to carry out a virtual ablation procedure and returns the
% ablation field as a binary array (1=ablated tissue) of length
% nnz(VoxelMat). VoxelMat is the cardiac geometry. int_scalar_mat is the
% scalar field visualized on the geometry during the session. This could be
% for example a lat map or a psi map. If the scalar field is constrained
% between 0 and 1 the function uses adas colormap. Otherwise gist_rainbow
% colormap (CARTO) is used with 8 levels (isochrones).
% res is the resolution and transmural the transmural coordinate, which can
% be passed as a one dimension array of length nnz(VoxelMat).
% The user can also provide an adas directory and adas filename to
% visualize adas corridors superimposed to the geometry. If these are nit
% provided the corridors visibility parametr in the interface is ignored.
% 
% 
[Ny,Nx,Nz]=size(VoxelMat);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
ind_sel=VoxelMat;
grx=gridx(ind_sel);
gry=gridy(ind_sel);
grz=gridz(ind_sel);

R=5;
R2=R^2;


[FV,extInd]=computeSurface(VoxelMat,res);
FV_cut=FV;
extInd_cut=extInd;
hf=figure;
hf.MenuBar='figure';
load('colormaps.mat','adas','gist_rainbow')
if max(int_scalar_mat)<=1
h=PlotVoxel(FV_cut,int_scalar_mat,extInd_cut,adas(1000));
else
    h=PlotVoxel(FV_cut,int_scalar_mat,extInd_cut,gist_rainbow(8));
end
ax=h.Parent;
hold on
handle_scatter=scatter3([],[],[],'filled','o','MarkerFaceColor','red');
handle_scatter_int=scatter3([],[],[],'filled','o','MarkerFaceColor','green');
ablation=false(nnz(VoxelMat),1);
plot_ablation=false(size(VoxelMat));
[FV1,extInd1]=computeSurface(plot_ablation,res);
FV_cut1=FV1;
extInd_cut1=extInd1;
h1=PlotVoxel(FV_cut1,extInd_cut1,plot_ablation,'y',1,0);
% h1.AlphaData=plot_ablation(extInd_cut);
h.ButtonDownFcn = {@click,handle_scatter,false};
if numel(transmural)==numel(VoxelMat)
    transmural=transmural(VoxelMat(:));
end
layer=9;
layer_thr=(layer-1)/8;
layer_select=false(size(VoxelMat));
layer_select(VoxelMat)=transmural<=layer_thr;
if nargin>5
corr=cell(9,1);
Points=[];
for ii=1:9
    VTK_shell=read_vtk([adas_folder '/' adas_filename '_LV_' num2str(ii*10) '_ENH_DE-MRI 2D.vtk']);
    Points=[Points; VTK_shell.points];
end

max_s=max(Points);
min_s=min(Points);
max_v=[max(grx) max(gry) max(grz)];
min_v=[min(grx) min(gry) min(grz)];

trasl=(min_s-min_v+max_s-max_v)/2;

for ii=1:9
    VTK_shell=read_vtk([adas_folder '/' adas_filename '_LV_' num2str(ii*10) '_CORR_DE-MRI 2D.vtk']);
    corr{ii}=patch('Faces',VTK_shell.elements,'Vertices',VTK_shell.points-trasl,'FaceColor','w','EdgeColor','none');
end

end

corr_v=1;

while true
    clc
    prev_abl=ablation;
    if corr_v==1
        corr_vis='all corridors';
    elseif corr_v==2
        corr_vis='only corridors on the current layer';
    elseif corr_v==3
        corr_vis='none';
    end
    str=input(['press:\n Q to exit \n'...
        ' L to add light and update \n' ...
        ' R to change the ablation sphere size (current R value: ' num2str(R) ')\n' ...
        ' B to delete last selected point \n'...
        ' C to define a cutplane \n'...
        ' F to flip the current cutplane \n'...
        ' Z to remove the current cutplane \n'...
        ' S to select ADAS layer (current: ' num2str(layer) ') \n' ...
        ' M to change color limits \n' ...
        ' V to change corridors visibility (current mode: ' corr_vis ') \n'...
        ' any other to update\n'],'s');
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

    if str=='B'
        if ~isempty(handle_scatter.XData)
            handle_scatter.XData=handle_scatter.XData(1:end-1);
            handle_scatter.YData=handle_scatter.YData(1:end-1);
            handle_scatter.ZData=[handle_scatter.ZData(1:end-1)];
        end
        continue
    end

    if str=='M'
        clc
        strLim=input('Insert color limits inside square brackets \n');
        ax.CLim=strLim;
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
        [Ny,Nx,Nz]=size(VoxelMat);
        x=(0.5:1:(Nx-0.5))*res;
        y=(0.5:1:(Ny-0.5))*res;
        z=(0.5:1:(Nz-0.5))*res;
        [gridx, gridy, gridz]=meshgrid(x,y,z);
        plane=normal(1)*gridx+normal(2)*gridy+normal(3)*gridz<d;
        [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane & layer_select,res);
        [FV_cut1,~]=computeSurface(plot_ablation & plane & layer_select,res);
        h.Faces=FV_cut.faces;
        h.Vertices=FV_cut.vertices;
        h.CData=int_scalar_mat(extInd_cut);
        h1.Faces=FV_cut1.faces;
        h1.Vertices=FV_cut1.vertices;
        handle_scatter_int.XData=[];
        handle_scatter_int.YData=[];
        handle_scatter_int.ZData=[];
        h.ButtonDownFcn = {@click,handle_scatter,false};
        continue
    end

    if str=='F'
        if exist('plane','var')
            plane=not(plane);
            [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane & layer_select,res);
            [FV_cut1,~]=computeSurface(plot_ablation & plane & layer_select,res);
            h.Faces=FV_cut.faces;
            h.Vertices=FV_cut.vertices;
            h.CData=int_scalar_mat(extInd_cut);
            h1.Faces=FV_cut1.faces;
            h1.Vertices=FV_cut1.vertices;
        end
        continue
    end

    if str=='Z'
        if exist('plane','var')
            [FV_cut,extInd_cut]=computeSurface(VoxelMat & layer_select,res);
            [FV_cut1,~]=computeSurface(plot_ablation & layer_select,res);
            clear plane
            h.Faces=FV_cut.faces;
            h.Vertices=FV_cut.vertices;
            h.CData=int_scalar_mat(extInd_cut);
            h1.Faces=FV_cut1.faces;
            h1.Vertices=FV_cut1.vertices;
        end
        continue
    end

    % if str=='M'
    %     clc
    %     strLim=input('Insert color limits inside square brackets \n');
    %     ax.CLim=strLim;
    %     continue
    % end

    if str=='S'
        layer=round(input('Select layer to view (from 1 to 9) \n'));
        layer_thr=(layer-1)/8;
        if layer_thr>1 || layer_thr<0
            warning('Layer should be an integer from 1 to 9')
            continue
        end
        layer_select=false(size(VoxelMat));
        layer_select(VoxelMat)=transmural<=layer_thr;
        if exist('plane','var')
            [FV_cut,extInd_cut]=computeSurface(VoxelMat & plane & layer_select,res);
            [FV_cut1,~]=computeSurface(plot_ablation & plane & layer_select,res);
        else
            [FV_cut,extInd_cut]=computeSurface(VoxelMat & layer_select,res);
            [FV_cut1,~]=computeSurface(plot_ablation & layer_select,res);
        end
        h.Faces=FV_cut.faces;
        h.Vertices=FV_cut.vertices;
        h.CData=int_scalar_mat(extInd_cut);
        h1.Faces=FV_cut1.faces;
        h1.Vertices=FV_cut1.vertices;
        if corr_v==2
            for ii=1:9
                try
                    corr{ii}.Visible=0;
                catch
                end
            end
            try
                corr{layer}.Visible=1;
            catch
            end
        end
        continue
    end

    if str=='V'
        str_v=input('Input visibility of corridors (1=all,2=current layer,3=none) \n');
        if str_v==1
            for ii=1:9
                try
                    corr{ii}.Visible=1;
                catch
                end
            end
            corr_v=str_v;
        elseif str_v==2
            for ii=1:9
                try
                    corr{ii}.Visible=0;
                catch
                end
            end
            try
                corr{layer}.Visible=1;
            catch
            end
            corr_v=str_v;
        elseif str_v==3
            for ii=1:9
                try
                    corr{ii}.Visible=0;
                catch
                end
            end
            corr_v=str_v;
        else
        end      
    end


    ntips=length(handle_scatter.XData);
    for i=1:ntips
        new_sel=(grx-handle_scatter.XData(i)).^2+(gry-handle_scatter.YData(i)).^2+(grz-handle_scatter.ZData(i)).^2<R2;
        ablation=or(ablation,new_sel);
    end
    plot_ablation(ind_sel)=ablation;
    handle_scatter.XData=[];
    handle_scatter.YData=[];
    handle_scatter.ZData=[];
    h.CData=int_scalar_mat(extInd_cut);
    if exist('plane','var')
        [FV_cut1,~]=computeSurface(plot_ablation & plane & layer_select,res);
    else
        [FV_cut1,~]=computeSurface(plot_ablation & layer_select,res);
    end
    h1.Faces=FV_cut1.faces;
    h1.Vertices=FV_cut1.vertices;
    if ntips>0
        str2=input('Confirm volume (D to delete)?','s');
        if str2=='D'
            ablation=prev_abl;
            plot_ablation(ind_sel)=ablation;
            h.CData=int_scalar_mat(extInd_cut);
            if exist('plane','var')
                [FV_cut1,~]=computeSurface(plot_ablation & plane,res);
            else
                [FV_cut1,~]=computeSurface(plot_ablation,res);
            end
            h1.Faces=FV_cut1.faces;
            h1.Vertices=FV_cut1.vertices;
        end
    end

end


end