function tag=tagAtria(VoxelMat,res,suppress_plot)
% tagAtria perform semiautomatic tag procedure of an atrial geometry.
% The user is asked to manually select the valve and vein ring openings with
% simple mouse clicking.
% tag=tagSurface(VoxelMat,res) returns an array tag of
% dimensions (numel(VoxelMat),1) containing the tag for each voxel. tag=0
% means the voxel is outside the geometry. Below it is the menaning of
% other tags:
% -1-->tricuspid valve ring
% -2-->inferior caval vein ring
% -3-->superior caval vein ring
% -4-->left pulmonary vein rings
% -5-->right pulmonary vein rings
% -6-->mitral valve ring
% -7-->right atrium endocardium
% -8-->left atrium endocardium
% -9-->epicardium
%The voxelized geometry (VoxelMat) is the only mandatory input. Default
% resolution is 0.25 mm.
% tag=tagAtria(VoxelMat,res,true) avoids showing the tag figure

if nargin<2
    res=0.25;
end

if nargin<3
    suppress_plot=false;
end

[FV,extInd]=computeSurface(VoxelMat,res);

ind=1:numel(VoxelMat);
ind_surf=unique(extInd);
isSurf=false(size(ind));
isSurf(ind_surf)=1;
tag=zeros(size(ind));

%tricuspid valve ring->1
disp('Select tricuspid valve ring...')
TriValve=selectSurface(VoxelMat,FV,extInd,res,2);
tag(TriValve)=1;

%inferior caval vein ring->2
disp('Select inferior caval vein ring...')
icv=selectSurface(VoxelMat,FV,extInd,res,1.5);
tag(icv)=2;

%superior caval vein ring->3
disp('Select superior caval vein ring...')
scv=selectSurface(VoxelMat,FV,extInd,res,1.5);
tag(scv)=3;

%left pulmonary vein rings->4
disp('Select left pulmonary vein rings...')
lpv=selectSurface(VoxelMat,FV,extInd,res,1.5);
tag(lpv)=4;

%right pulmonary vein rings->5
disp('Select right pulmonary vein rings...')
rpv=selectSurface(VoxelMat,FV,extInd,res,1.5);
tag(rpv)=5;

%mitral valve->6
disp('Select mitral valve ring...')
mv=selectSurface(VoxelMat,FV,extInd,res,2);
tag(mv)=6;

disp('Done!')

%% automatic distinguish right/left endo and epi surfaces
[Ny,Nx,~]=size(VoxelMat);
tagSurf=tag(ind_surf);
indToTag=find(tagSurf==0);

kk=randi(length(indToTag));
new_seed=ind(ind_surf(indToTag(kk)));
tag(new_seed)=7;

Nsurf=length(ind_surf);

seed=new_seed;
while nnz(tag)<Nsurf

    nn=length(seed);

    new_seed=[];



    for jj=1:nn


        %+dy
        fw_y=seed(jj)+1;
        if isSurf(fw_y)
            if tag(fw_y)==0
                tag(fw_y)=tag(seed(jj));
                new_seed=[new_seed; fw_y];
            end
        end

        %-dy
        bw_y=seed(jj)-1;
        if isSurf(bw_y)
            if tag(bw_y)==0

                tag(bw_y)=tag(seed(jj));
                new_seed=[new_seed; bw_y];

            end
        end

        %+dx
        fw_x=seed(jj)+Ny;
        if isSurf(fw_x)
            if tag(fw_x)==0
                tag(fw_x)=tag(seed(jj));
                new_seed=[new_seed; fw_x];

            end
        end

        %-dx
        bw_x=seed(jj)-Ny;
        if isSurf(bw_x)
            if tag(bw_x)==0

                tag(bw_x)=tag(seed(jj));
                new_seed=[new_seed; bw_x];

            end
        end

        %+dz
        fw_z=seed(jj)+Ny*Nx;
        if isSurf(fw_z)
            if tag(fw_z)==0

                tag(fw_z)=tag(seed(jj));
                new_seed=[new_seed; fw_z];

            end
        end


        %-dz
        bw_z=seed(jj)-Ny*Nx;
        if isSurf(bw_z)
            if tag(bw_z)==0

                tag(bw_z)=tag(seed(jj));
                new_seed=[new_seed; bw_z];

            end
        end


        %+dx +dy
        fx_fy=seed(jj)+Ny+1;
        if isSurf(fx_fy)
            if tag(fx_fy)==0

                tag(fx_fy)=tag(seed(jj));
                new_seed=[new_seed; fx_fy];

            end
        end

        %-dx -dy
        bx_by=seed(jj)-Ny-1;
        if isSurf(bx_by)
            if tag(bx_by)==0


                tag(bx_by)=tag(seed(jj));
                new_seed=[new_seed; bx_by];

            end
        end

        %-dx +dy
        bx_fy=seed(jj)-Ny+1;
        if isSurf(bx_fy)
            if tag(bx_fy)==0

                tag(bx_fy)=tag(seed(jj));
                new_seed=[new_seed; bx_fy];

            end
        end



        %+dx -dy
        fx_by=seed(jj)+Ny-1;
        if isSurf(fx_by)
            if tag(fx_by)==0

                tag(fx_by)=tag(seed(jj));
                new_seed=[new_seed; fx_by];

            end
        end


        %+dx +dz
        fx_fz=seed(jj)+Ny+Ny*Nx;
        if isSurf(fx_fz)
            if tag(fx_fz)==0

                tag(fx_fz)=tag(seed(jj));
                new_seed=[new_seed; fx_fz];
            end
        end

        %-dx -dz
        bx_bz=seed(jj)-Ny-Ny*Nx;
        if isSurf(bx_bz)
            if tag(bx_bz)==0

                tag(bx_bz)=tag(seed(jj));
                new_seed=[new_seed; bx_bz];

            end
        end

        %+dx -dz
        fx_bz=seed(jj)+Ny-Ny*Nx;
        if isSurf(fx_bz)
            if tag(fx_bz)==0

                tag(fx_bz)=tag(seed(jj));
                new_seed=[new_seed; fx_bz];

            end
        end

        %-dx +dz
        bx_fz=seed(jj)-Ny+Ny*Nx;
        if isSurf(bx_fz)
            if tag(bx_fz)==0

                tag(bx_fz)=tag(seed(jj));
                new_seed=[new_seed; bx_fz];

            end
        end

        %+dy %+dz
        fy_fz=seed(jj)+1+Ny*Nx;
        if isSurf(fy_fz)
            if tag(fy_fz)==0

                tag(fy_fz)=tag(seed(jj));
                new_seed=[new_seed; fy_fz];

            end
        end

        %+dy -dz
        fy_bz=seed(jj)+1-Ny*Nx;
        if isSurf(fy_bz)
            if tag(fy_bz)==0

                tag(fy_bz)=tag(seed(jj));
                new_seed=[new_seed; fy_bz];

            end
        end


        %-dy +dz
        by_fz=seed(jj)-1+Ny*Nx;
        if isSurf(by_fz)
            if tag(by_fz)==0

                tag(by_fz)=tag(seed(jj));
                new_seed=[new_seed; by_fz];

            end
        end

        %-dy -dz
        by_bz=seed(jj)-1-Ny*Nx;
        if isSurf(by_bz)
            if tag(by_bz)==0

                tag(by_bz)=tag(seed(jj));
                new_seed=[new_seed; by_bz];

            end
        end



    end




    if isempty (new_seed)

        tagSurf=tag(ind_surf);
        indToTag=find(tagSurf==0);
        kk=randi(length(indToTag));
        new_seed=ind(ind_surf(indToTag(kk)));
        tag(new_seed)=tag(seed(end))+1;
    end


    seed=new_seed;

end
ntags=max(tag);
if ntags<7
    warning('User-selected rings does not completely separate endocardial and epicardial surfaces. Consider to run again the function to newly select the base')
else
    nel_tag=zeros(ntags,1);
    for i=7:ntags
        nel_tag(i)=nnz(tag==i);
    end
    [~,sorted_tag]=sort(nel_tag);
    tmp=tag;
    tmp(tag==sorted_tag(end))=9; %epicardium
    tmp(tag==sorted_tag(end-2))=8; %endocardium LA
    tmp(tag==sorted_tag(end-1))=7; %endocardium RA
    tag=tmp;
    toExc=find(tag>9);
    while ~isempty(toExc)
        for ii=[1 -1 Ny -Ny Ny*Nx -Ny*Nx 1+Ny -1-Ny 1-Ny -1+Ny Ny+Ny*Nx -Ny-Nx*Ny Ny-Ny*Nx -Ny+Ny*Nx 1+Ny*Nx -1-Ny*Nx 1-Ny*Nx -1+Ny*Nx]
            upt=tag(toExc+ii)<=9 & tag(toExc+ii)>0;
            tag(toExc(upt))=tag(toExc(upt)+ii);
            toExc(upt)=[];
            if isempty(toExc)
                break;
            end
        end
    end
end
tag=tag';

if not(suppress_plot)
    figure;
    PlotVoxel(FV,tag,extInd);
end
