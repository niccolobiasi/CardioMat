function tag=tagSurface(VoxelMat,res,suppress_plot)
% tagSurface selects the basal, epicardial, and endocardial surfaces in a
% left ventricle based on user inputs. The user is asked to manually select
% the base with simple mouse clicking.
% tag=tagSurface(VoxelMat,res) returns an array tag of
% dimensions (numel(VoxelMat),1) containing the tag for each voxel. tag=0
% means the voxel is outside the geometry. tag=1 indicates the
% user-selected base. tag>1 are used for automatically detected surfaces.
% tag=4 indicates the epicardial surface. tag=3 indicates LV endocardium.
% tag=2 indicates RV endocardium.
% If monoventricular mesh is provided, the function assumes it is LV
% geometry with septum.
% The tag are assigned based on the surface area: the largest is
% assigned to the epicardium and the others to the endocardium.
% The largest between the two endocardial surface is considered RV.
% If the user-selected base does not completely separate endo and epi surfaces a
% warning is raised and all the surface is considered endocardial.
% The voxelized geometry (VoxelMat) is the only mandatory input. Default
% resolution is 0.25 mm.
% tag=tagSurface(VoxelMat,res,true) avoids showing the tag figure.



if nargin<2
    res=0.25;
end

if nargin<3
    suppress_plot=false;
end

[Ny,Nx,~]=size(VoxelMat);
ind=1:numel(VoxelMat);

[FV,extInd]=computeSurface(VoxelMat,res);
ind_surf=unique(extInd);
isSurf=false(size(ind));
isSurf(ind_surf)=1;

disp('Select base surface')
plot_selected=selectSurface(VoxelMat,FV,extInd,res,5);
%%
tag=zeros(size(ind));
tag(plot_selected)=1;
tagSurf=tag(ind_surf);
indToTag=find(tagSurf==0);

kk=randi(length(indToTag));
new_seed=ind(ind_surf(indToTag(kk)));
tag(new_seed)=2;

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
        tag(new_seed)=tag(seed(1))+1;
    end


    seed=new_seed;


end

%%
ntags=max(tag);
if ntags<3
    warning('User-selected base does not completely separate endocardial and epicardial surfaces. Consider to run again the function to newly select the base')
else
    nel_tag=zeros(ntags,1);
    for i=2:ntags
        nel_tag(i)=nnz(tag==i);
    end
    [~,sorted_tag]=sort(nel_tag);
    tmp=double(isSurf(:));
    tmp(tag==sorted_tag(end))=4; %epicardium
    if ntags>3
        if nnz(tag==sorted_tag(end-2))>nnz(plot_selected)
            tmp(tag==sorted_tag(end-2))=3; %endocardium LV
            tmp(tag==sorted_tag(end-1))=2; %endocardium RV
        else
            tmp(tag==sorted_tag(end-1))=3; %endocardium LV
        end
    else
        tmp(tag==sorted_tag(end-1))=3; %endocardium LV
    end
    tag=tmp;
end


if not(suppress_plot)
    figure;
    PlotVoxel(FV,tag,extInd);
end


end