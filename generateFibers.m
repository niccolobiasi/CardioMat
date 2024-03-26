function [f,s,t]=generateFibers(VoxelMat,apicobasal,transmural,options)
% generateFibers generates longitudinal (f), normal-to-sheet (s), and
% transverse (t) fiber directions for a voxelized cardiac geometry
% VoxelMat by using the method described by Bayer et al.. Transmural and apicobasal coordinates are reuquired as input.
% An additional input structure options is used to pass parameter values.
% options may contain the following fields:
%- res: specify the resolution (used only for plot, default=0.25 mm)
%- plot: specify if create a plot or not (default= false)
%- alpha_epi_lv: specify the alpha angle at the epicardium in the left 
% ventricle (default=-60)
%- alpha_endo_lv: specify the alpha angle at the endocardium in the left
% ventricle (default=60)
%- alpha_epi_rv: specify the alpha angle at the epicardium in the right 
% ventricle (default=-25)
%- alpha_endo_rv: specify the alpha angle at the endocardium in the right
% ventricle (default=90)
%- beta_epi_lv: specify the beta angle at the epicardium in the left
%ventricle (default=20)
%- beta_endo_lv: specify the beta angle at the endocardium in the left
%ventricle (default=-20)
%- beta_epi_rv: specify the beta angle at the epicardium in the right
%ventricle (default=20)
%- beta_endo_rv: specify the beta angle at the endocardium in the right
%ventricle (default=0)
%- transmural_bv: a transmural field for biventricular geometries (pointing
% from left to right ventricle in the septum).
%
% f=generateFibers(VoxelMat,apicobasal,transmural,options) returns a
% structure f with field 'x','y', and 'z' containing the components of the
% longitudinal fiber direction in the three cartesian directions for each 
% voxel belonging to the voxelized geometry.
%
% [f,s,t]=generateFibers(VoxelMat,apicobasal,transmural,options) also
% returns normal-to-sheet and transverse fiber directions.

if nargin<4
    options=[];
end

if isfield(options,'res')
    res=options.res;
else
    res=0.25;
end

if isfield(options,'plot')
    doplot=options.plot;
else
    doplot=false;
end

if isfield(options,'aplha_epi_lv')
    alpha_epi_lv=options.alpha_epi_lv;
else 
    alpha_epi_lv=deg2rad(-60);
end

if isfield(options, 'alpha_endo_lv')
    alpha_endo_lv=options.alpha_endo_lv;
else
    alpha_endo_lv=deg2rad(60);
end

if isfield(options,'aplha_epi_rv')
    alpha_epi_rv=options.alpha_epi_rv;
else 
    alpha_epi_rv=deg2rad(-25);
end

if isfield(options, 'alpha_endo_rv')
    alpha_endo_rv=options.alpha_endo_rv;
else
    alpha_endo_rv=deg2rad(90);
end

if isfield(options, 'beta_epi_lv')
    beta_epi_lv=options.beta_epi_lv;
else
    beta_epi_lv=deg2rad(20);
end

if isfield(options,'beta_endo_lv')
    beta_endo_lv=options.beta_endo_lv;
else
    beta_endo_lv=deg2rad(-20);
end

if isfield(options, 'beta_epi_rv')
    beta_epi_rv=options.beta_epi_rv;
else
    beta_epi_rv=deg2rad(20);
end

if isfield(options,'beta_endo_rv')
    beta_endo_rv=options.beta_endo_rv;
else
    beta_endo_rv=deg2rad(0);
end

if isfield(options,'transmural_bv')
    transmural_fibers=options.transmural_bv;
else
    transmural_fibers=transmural;
end

if isfield(options,'interventricular')
    interventricular=options.interventricular;
else
    interventricular=zeros(numel(VoxelMat),1);
end

S=gradMat(VoxelMat,'x');
tmx=S*transmural_fibers(VoxelMat);
abx=S*apicobasal(VoxelMat);
S=gradMat(VoxelMat,'y');
tmy=S*transmural_fibers(VoxelMat);
aby=S*apicobasal(VoxelMat);
S=gradMat(VoxelMat,'z');
tmz=S*transmural_fibers(VoxelMat);
abz=S*apicobasal(VoxelMat);

%%
norm_ab=sqrt(abx.^2+aby.^2+abz.^2);
elx=abx./(norm_ab);
ely=aby./(norm_ab);
elz=abz./(norm_ab);

tmp=elx.*tmx+ely.*tmy+elz.*tmz;
etx=tmx-tmp.*elx;
ety=tmy-tmp.*ely;
etz=tmz-tmp.*elz;
norm_et=sqrt(etx.^2+ety.^2+etz.^2);
etx=etx./(norm_et);
ety=ety./(norm_et);
etz=etz./(norm_et);

ecx=ely.*etz-elz.*ety;
ecy=elz.*etx-elx.*etz;
ecz=elx.*ety-ely.*etx;
norm_ec=sqrt(ecx.^2+ecy.^2+ecz.^2);
ecx=ecx./(norm_ec);
ecy=ecy./(norm_ec);
ecz=ecz./(norm_ec);
%%

elx_plot=zeros(size(VoxelMat));
elx_plot(VoxelMat)=elx;
ely_plot=zeros(size(VoxelMat));
ely_plot(VoxelMat)=ely;
elz_plot=zeros(size(VoxelMat));
elz_plot(VoxelMat)=elz;

etx_plot=zeros(size(VoxelMat));
etx_plot(VoxelMat)=etx;
ety_plot=zeros(size(VoxelMat));
ety_plot(VoxelMat)=ety;
etz_plot=zeros(size(VoxelMat));
etz_plot(VoxelMat)=etz;

ecx_plot=zeros(size(VoxelMat));
ecx_plot(VoxelMat)=ecx;
ecy_plot=zeros(size(VoxelMat));
ecy_plot(VoxelMat)=ecy;
ecz_plot=zeros(size(VoxelMat));
ecz_plot(VoxelMat)=ecz;
%%

% correct NaN
bad=isnan(elx_plot) | isnan(etx_plot) | isnan(ecx_plot);
ind_bad=find(bad);
[Ny,Nx,~]=size(VoxelMat);
ind_in=find(VoxelMat);
while not(isempty(ind_bad))
    for ii=[1 -1 Ny -Ny Ny*Nx -Ny*Nx]
        upt_bad=not(bad(ind_bad+ii)) & VoxelMat(ind_bad+ii);
        upt=ind_bad(upt_bad);

        elx_plot(upt)=elx_plot(upt+ii);
        ely_plot(upt)=ely_plot(upt+ii);
        elz_plot(upt)=elz_plot(upt+ii);

        etx_plot(upt)=etx_plot(upt+ii);
        ety_plot(upt)=ety_plot(upt+ii);
        etz_plot(upt)=etz_plot(upt+ii);

        ecx_plot(upt)=ecx_plot(upt+ii);
        ecy_plot(upt)=ecy_plot(upt+ii);
        ecz_plot(upt)=ecz_plot(upt+ii);

        bad=isnan(elx_plot) | isnan(etx_plot) | isnan(ecx_plot);
        ind_bad=find(bad);
    end

end
elx=elx_plot(ind_in);
ely=ely_plot(ind_in);
elz=elz_plot(ind_in);

etx=etx_plot(ind_in);
ety=ety_plot(ind_in);
etz=etz_plot(ind_in);

ecx=ecx_plot(ind_in);
ecy=ecy_plot(ind_in);
ecz=ecz_plot(ind_in);
%%

alpha=(alpha_endo_lv*(1-transmural(ind_in))+alpha_epi_lv*transmural(ind_in)).*(1-interventricular(ind_in))+...
    (alpha_endo_rv*(1-transmural(ind_in))+alpha_epi_rv*transmural(ind_in)).*interventricular(ind_in);


beta=(beta_endo_lv*(1-transmural(ind_in))+beta_epi_lv*transmural(ind_in)).*(1-interventricular(ind_in))+...
    (beta_endo_rv*(1-transmural(ind_in))+beta_epi_rv*transmural(ind_in)).*interventricular(ind_in);


cos_alpha=cos(alpha);
sin_alpha=sin(alpha);

cos_beta=cos(beta);
sin_beta=sin(beta);

f.x=ecx.*cos_alpha+elx.*sin_alpha;
f.y=ecy.*cos_alpha+ely.*sin_alpha;
f.z=ecz.*cos_alpha+elz.*sin_alpha;

if nargout>1
    efx=-ecx.*sin_alpha+elx.*cos_alpha;
    efy=-ecy.*sin_alpha+ely.*cos_alpha;
    efz=-ecz.*sin_alpha+elz.*cos_alpha;

    s.x=efx.*cos_beta-etx.*sin_beta;
    s.y=efy.*cos_beta-ety.*sin_beta;
    s.z=efz.*cos_beta-etz.*sin_beta;
end

if nargout>2
    t.x=efx.*sin_beta+etx.*cos_beta;
    t.y=efy.*sin_beta+ety.*cos_beta;
    t.z=efz.*sin_beta+etz.*cos_beta;
end

if doplot==1
    options_plot.res=res;
    plotFibers(f.x,f.y,f.z,VoxelMat,transmural,options_plot);
end

end
