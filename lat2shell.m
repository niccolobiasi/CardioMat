function lat2shell(folder_out,adas_folder,adas_filename,filename_geom,act_times,show_scar,ablation,res,levels,color_name)
% lat2shell maps local activation time array act_times to the surface meshes
% generated by ADAS.
% lat2shell(folder_out,adas_folder,adas_filename,filename_geom,act_times)
% generates a folder_out folder with 9 vtk file named "lats_ii.vtk". Each
% vtk file corresponds to an ADAS shell. ADAS .vtk shells with filename
% adas_filename in the folder adas_folder should be in the format:
% "adas_filename"_LV_ii_ENH_DE-MRI 2D.vtk" where ii goes from 10 to 90.
% The user should also specify a mat file containing the geometry features
% of the model.
% The user can also specify if showing core scar tissue in white in the
% exported maps. Additionally, a further scar tissue mask can be passed as
% a vector of length nnz(VoxelMat).
% The user can also specify the resolution of the model( default: 0.25mm),
% the number of color levels to be used (default:256), and the name of the
% colormap to choose between default matlab colormaps and additional
% colormaps in colormaps.mat file (default: gist_rainbow).

if iscell(act_times)
    act_times=extractLat(act_times,'first');
end

load([filename_geom '.mat'],"ind_in",'VoxelMat');

if nargin<10
    color_name='gist_rainbow';
elseif isempty(color_name)
    color_name='gist_rainbow';
end

if nargin<9
    levels=256;
elseif isempty(levels)
    levels=256;
end

if nargin<8
    res=0.25;
elseif isempty(res)
    res=0.25;
end

if exist(color_name,'file')
    cmap=feval(color_name,levels);
else
    load('colormaps.mat',color_name);
    if ~exist(color_name,"var")
        warning('using default colormap')
        color_name='gist_rainbow';
        load('colormaps.mat',color_name);
    end
    cmap_h=eval(color_name);
    cmap=cmap_h(levels);
end

if nargin<7
    ablation=false(nnz(VoxelMat),1);
elseif isempty(ablation)
    ablation=false(nnz(VoxelMat),1);
end

if nargin<6
    show_scar=false;
elseif isempty(show_scar)
    show_scar=false;
end


[Ny,Nx,Nz]=size(ind_in);
range_x=(Nx-1)*res;
range_y=(Ny-1)*res;
range_z=(Nz-1)*res;


VTK_shells=cell(9,1);
Points=[];
Scalars=[];

for ii=1:9
    VTK_shells{ii}=read_vtk([adas_folder '/' adas_filename '_LV_' num2str(ii*10) '_ENH_DE-MRI 2D.vtk']);
    Points=[Points; VTK_shells{ii}.points];
    Scalars=[Scalars; VTK_shells{ii}.scalar];
end

Xmax=max(Points(:,1));
Ymax=max(Points(:,2));
Zmax=max(Points(:,3));

Xmin=min(Points(:,1));
Ymin=min(Points(:,2));
Zmin=min(Points(:,3));

bound_x=(range_x-(Xmax-Xmin))/2;
bound_y=(range_y-(Ymax-Ymin))/2;
bound_z=(range_z-(Zmax-Zmin))/2;
x=(Xmin-bound_x):res:(Xmax+bound_x);
y=(Ymin-bound_y):res:(Ymax+bound_y);
z=(Zmin-bound_z):res:(Zmax+bound_z);

[gridx, gridy,gridz]=meshgrid(x,y,z);

tokeep= not(isnan(act_times) | isinf(act_times));
good=ind_in;
good(ind_in)=tokeep;
[~,extInd]=computeSurface(VoxelMat,res);
toFill=~isnan(extrapField(ones(nnz(VoxelMat),1),VoxelMat,extInd,[],30));
toFill=reshape(toFill,size(VoxelMat));
[~,extInd]=computeSurface(good,res);
act_times=extrapField(act_times(tokeep),good,extInd,toFill);
act_times=reshape(act_times,size(ind_in));
if show_scar
    sth=2;
else
    sth=5;
end
act_times=imgaussfilt3(act_times,sth);
P = [2 1 3];
X = permute(gridx, P);
Y = permute(gridy, P);
Z = permute(gridz, P);
V = permute(act_times, P);
F = griddedInterpolant(X,Y,Z,V);

if show_scar
    [~,extInd]=computeSurface(VoxelMat,res);
    fib=extrapField(ablation,VoxelMat,extInd,[],30);
    fib=reshape(fib,size(VoxelMat));
    V_fib=permute(fib,P);
    F_fib=griddedInterpolant(X,Y,Z,V_fib,'nearest','nearest');
end

lats=F(Points);
if show_scar
    scar_all=logical(F_fib(Points));
end
mkdir(folder_out)
curr=1;
if show_scar
    cmap=[cmap; 1 1 1];
end
for ii=1:9
    nn=length(VTK_shells{ii}.points);
    if show_scar
        scar=VTK_shells{ii}.scalar>(1-1e-6) | scar_all(curr:curr+nn-1);
    else
        scar=false(nn,1);
    end
    shell_lat=lats(curr:curr+nn-1);
    norm_lat=levels/(levels+1)*(shell_lat-min(shell_lat(~scar)))/(max(shell_lat(~scar))-min(shell_lat(~scar)))-1e-3;
    norm_lat(scar)=1;
    write_vtk_polydata([folder_out '/lats_' num2str(ii*10) '.vtk'], ...
        VTK_shells{ii}.points,VTK_shells{ii}.elements, norm_lat,cmap)

    curr=curr+nn;
end