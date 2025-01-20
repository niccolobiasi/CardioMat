function fields_RA=computeFieldsRA(RA,tag,res)
% fields_RA=computeFieldsRA(RA,tag,res) returns a cell array of 5 elements containing the
% distance fields needed for the assignement in the fiber orientation in
% the right atrium (voxelized geometry RA). 
% For further information see the related paper by Piersanti et al. (2021)
% The function requires the surface tags tag as produced by tagAtria
% function. The resolution is optional (default: 0.25 mm). Manual selection
% of the right atrial appendage is required.

if nargin<3 
    res=0.25;
end


[FV_r,extInd_r]=computeSurface(RA,res);
ind_surf_r=unique(extInd_r);
isSurf_r=false(size(RA));
isSurf_r(ind_surf_r)=1;
eik=solveEikonal(RA,find((tag==9 | (tag==0 & isSurf_r(:))) & RA(:)));
eik1=solveEikonal(RA,find(tag==7));
phi_r=eik1./(eik+eik1);
eik_icv=solveEikonal(RA,find(tag==2));
eik_tv=solveEikonal(RA,find(tag==1));
eik_scv=solveEikonal(RA,find(tag==3));
disp('Select right atrial appendage...')
ra_app=selectSurface(RA,FV_r,extInd_r,res,1, true);
disp('Done!')
eik_rapp=solveEikonal(RA,find(ra_app));

eik_tv_scv_rapp=solveEikonal(RA,find(tag==1| tag==3 | ra_app(:)));
eik_icv_scv_rapp=solveEikonal(RA, find(tag==2 | tag==3 |ra_app(:)));
eik_icv_tv_scv=solveEikonal(RA, find(tag==1| tag==2| tag==3));
eik_icv_tv_rapp=solveEikonal(RA, find(tag==2| tag==1| ra_app(:)));

psi_ab_ra=(3*(eik_rapp./(eik_rapp+eik_icv_tv_scv))+1*(eik_scv./(eik_scv+eik_icv_tv_rapp))...
    -1*(eik_tv./(eik_tv+eik_icv_scv_rapp))-3*(eik_icv./(eik_icv+eik_tv_scv_rapp))+1)/2;

eik2=solveEikonal(RA, find(tag==3 | ra_app(:)));
psi_v_ra=eik2./(eik2+eik_icv);

indTV=find(tag==1);
[yy,xx,zz]=ind2sub(size(RA),indTV);
centerMV=[mean(xx) mean(yy) mean(zz)]*res;

indICV=find(tag==2);
[yy,xx,zz]=ind2sub(size(RA),indICV);
centerICV=[mean(xx) mean(yy) mean(zz)]*res;

indSCV=find(tag==3);
[yy,xx,zz]=ind2sub(size(RA),indSCV);
centerSCV=[mean(xx) mean(yy) mean(zz)]*res;


p1=centerMV;
p2=centerICV;
p3=centerSCV;

normal = cross(p1 - p2, p1 - p3);
d = (p1(1)*normal(1) + p1(2)*normal(2) + p1(3)*normal(3));

[Ny,Nx,Nz]=size(RA);
x=(0.5:1:(Nx-0.5))*res;
y=(0.5:1:(Ny-0.5))*res;
z=(0.5:1:(Nz-0.5))*res;
[gridx, gridy, gridz]=meshgrid(x,y,z);
gridx=gridx(RA(:));
gridy=gridy(RA(:));
gridz=gridz(RA(:));

vv=normal(1)*gridx+normal(2)*gridy+normal(3)*gridz-d;
planeINT=0.5>abs(vv)/norm(normal);
plane=false(Ny,Nx,Nz);
plane(RA(:))=planeINT;

plane_epi=reshape(plane(:) & (tag==9 | tag==2 | tag==3),Ny,Nx,Nz);

% automatic extraction of top band
seed=find(plane_epi(:) & tag==3);
[~,seed_m]=max(psi_ab_ra(seed));
seed=seed(seed_m);

top=false(Ny,Nx,Nz);

while ~isempty(seed)

    nn=length(seed);

    new_seed=[];



    for jj=1:nn


        %+dy
        fw_y=seed(jj)+1;
        if plane_epi(fw_y)
            if top(fw_y)==0
                top(fw_y)=1;
                new_seed=[new_seed; fw_y];
            end
        end

        %-dy
        bw_y=seed(jj)-1;
        if plane_epi(bw_y)
            if top(bw_y)==0

                top(bw_y)=1;
                new_seed=[new_seed; bw_y];

            end
        end

        %+dx
        fw_x=seed(jj)+Ny;
        if plane_epi(fw_x)
            if top(fw_x)==0
                top(fw_x)=1;
                new_seed=[new_seed; fw_x];

            end
        end

        %-dx
        bw_x=seed(jj)-Ny;
        if plane_epi(bw_x)
            if top(bw_x)==0

                top(bw_x)=1;
                new_seed=[new_seed; bw_x];

            end
        end

        %+dz
        fw_z=seed(jj)+Ny*Nx;
        if plane_epi(fw_z)
            if top(fw_z)==0

                top(fw_z)=1;
                new_seed=[new_seed; fw_z];

            end
        end


        %-dz
        bw_z=seed(jj)-Ny*Nx;
        if plane_epi(bw_z)
            if top(bw_z)==0

                top(bw_z)=1;
                new_seed=[new_seed; bw_z];

            end
        end


        %+dx +dy
        fx_fy=seed(jj)+Ny+1;
        if plane_epi(fx_fy)
            if top(fx_fy)==0

                top(fx_fy)=1;
                new_seed=[new_seed; fx_fy];

            end
        end

        %-dx -dy
        bx_by=seed(jj)-Ny-1;
        if plane_epi(bx_by)
            if top(bx_by)==0


                top(bx_by)=1;
                new_seed=[new_seed; bx_by];

            end
        end

        %-dx +dy
        bx_fy=seed(jj)-Ny+1;
        if plane_epi(bx_fy)
            if top(bx_fy)==0

                top(bx_fy)=1;
                new_seed=[new_seed; bx_fy];

            end
        end



        %+dx -dy
        fx_by=seed(jj)+Ny-1;
        if plane_epi(fx_by)
            if top(fx_by)==0

                top(fx_by)=1;
                new_seed=[new_seed; fx_by];

            end
        end


        %+dx +dz
        fx_fz=seed(jj)+Ny+Ny*Nx;
        if plane_epi(fx_fz)
            if top(fx_fz)==0

                top(fx_fz)=1;
                new_seed=[new_seed; fx_fz];
            end
        end

        %-dx -dz
        bx_bz=seed(jj)-Ny-Ny*Nx;
        if plane_epi(bx_bz)
            if top(bx_bz)==0

                top(bx_bz)=1;
                new_seed=[new_seed; bx_bz];

            end
        end

        %+dx -dz
        fx_bz=seed(jj)+Ny-Ny*Nx;
        if plane_epi(fx_bz)
            if top(fx_bz)==0

                top(fx_bz)=top(seed(jj));
                new_seed=[new_seed; fx_bz];

            end
        end

        %-dx +dz
        bx_fz=seed(jj)-Ny+Ny*Nx;
        if plane_epi(bx_fz)
            if top(bx_fz)==0

                top(bx_fz)=1;
                new_seed=[new_seed; bx_fz];

            end
        end

        %+dy %+dz
        fy_fz=seed(jj)+1+Ny*Nx;
        if plane_epi(fy_fz)
            if top(fy_fz)==0

                top(fy_fz)=1;
                new_seed=[new_seed; fy_fz];

            end
        end

        %+dy -dz
        fy_bz=seed(jj)+1-Ny*Nx;
        if plane_epi(fy_bz)
            if top(fy_bz)==0

                top(fy_bz)=1;
                new_seed=[new_seed; fy_bz];

            end
        end


        %-dy +dz
        by_fz=seed(jj)-1+Ny*Nx;
        if plane_epi(by_fz)
            if top(by_fz)==0

                top(by_fz)=1;
                new_seed=[new_seed; by_fz];

            end
        end

        %-dy -dz
        by_bz=seed(jj)-1-Ny*Nx;
        if plane_epi(by_bz)
            if top(by_bz)==0

                top(by_bz)=1;
                new_seed=[new_seed; by_bz];

            end
        end



    end


    seed=new_seed;


end
eik_top=solveEikonal(RA,find(top));
psi_r_ra=eik_top./(eik_tv+eik_top);

planeINT=vv>0;
plane=false(Ny,Nx,Nz);
plane(RA(:))=planeINT;
tv_septum=plane(:) & tag==1;
tv_free=~plane(:) & tag==1;

eik_tv_free=solveEikonal(RA,find(tv_free));
eik_tv_septum=solveEikonal(RA, find(tv_septum));
psi_w_ra=2*eik_tv_free./(eik_tv_free+eik_tv_septum)-1;
fields_RA={phi_r,psi_ab_ra,psi_v_ra, psi_r_ra, psi_w_ra};