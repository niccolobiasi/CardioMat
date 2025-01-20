function [f,s,n]=generateFiberLA(LA,fields_LA,par)
% generateFiberLA(LA,fields_LA,par) generate fiber orientation for the
% voxelized left-atrial geometry LA by using a rule-based method (see
% Piersanti et al. 2021). The function uses the distance fields fields_LA
% as computed by computeFieldsLA. Optionally, a par structure can be input
% defining the parameter values for bundles selection with subfields:
%   - tau_mv (default:0.85)
%   - tau_lpv (default: 0.85)
%   - tau_rpv (default: 0.2)


if nargin<3
    par=[];
end
if isfield(par,'tau_mv')
    tau_mv=par.tau_mv;
else
    tau_mv=0.85;
end

if isfield(par,'tau_lpv')
    tau_lpv=par.tau_lpv;
else
    tau_lpv=0.85;
end

if isfield(par,'tau_rpv')
    tau_rpv=par.tau_rpv;
else
    tau_rpv=0.2;
end

phi_l=fields_LA{1};
psi_ab_la=fields_LA{2};
psi_v_la=fields_LA{3};
psi_r_la=fields_LA{4};

S=gradMat(LA,'x');
tmx=S*phi_l(LA);
r_l_x=S*psi_r_la(LA);
v_l_x=S*psi_v_la(LA);
ab_l_x=S*psi_ab_la(LA);
S=gradMat(LA,'y');
tmy=S*phi_l(LA);
r_l_y=S*psi_r_la(LA);
v_l_y=S*psi_v_la(LA);
ab_l_y=S*psi_ab_la(LA);
S=gradMat(LA,'z');
tmz=S*phi_l(LA);
r_l_z=S*psi_r_la(LA);
v_l_z=S*psi_v_la(LA);
ab_l_z=S*psi_ab_la(LA);

norm_et=sqrt(tmx.^2+tmy.^2+tmz.^2);
etx=tmx./(norm_et);
ety=tmy./(norm_et);
etz=tmz./(norm_et);

kx=zeros(size(etx));
ky=zeros(size(etx));
kz=zeros(size(etx));
ind1=psi_r_la(LA)>=tau_mv;
kx(ind1)=r_l_x(ind1);
ky(ind1)=r_l_y(ind1);
kz(ind1)=r_l_z(ind1);

ind2=~ind1 & (psi_v_la(LA)>=tau_lpv | psi_v_la(LA)<=tau_rpv);
kx(ind2)=v_l_x(ind2);
ky(ind2)=v_l_y(ind2);
kz(ind2)=v_l_z(ind2);

ind3=~ind1 & ~ind2;
kx(ind3)=ab_l_x(ind3);
ky(ind3)=ab_l_y(ind3);
kz(ind3)=ab_l_z(ind3);

tmp=kx.*etx+ky.*ety+kx.*etz;
enx=kx-tmp.*etx;
eny=ky-tmp.*ety;
enz=kz-tmp.*etz;
norm_en=sqrt(enx.^2+eny.^2+enz.^2);
enx=enx./(norm_en);
eny=eny./(norm_en);
enz=enz./(norm_en);

elx=eny.*etz-enz.*ety;
ely=enz.*etx-enx.*etz;
elz=enx.*ety-eny.*etx;
norm_el=sqrt(elx.^2+ely.^2+elz.^2);
elx=elx./(norm_el);
ely=ely./(norm_el);
elz=elz./(norm_el);

elx_plot=zeros(size(LA));
elx_plot(LA)=elx;
ely_plot=zeros(size(LA));
ely_plot(LA)=ely;
elz_plot=zeros(size(LA));
elz_plot(LA)=elz;

etx_plot=zeros(size(LA));
etx_plot(LA)=etx;
ety_plot=zeros(size(LA));
ety_plot(LA)=ety;
etz_plot=zeros(size(LA));
etz_plot(LA)=etz;

enx_plot=zeros(size(LA));
enx_plot(LA)=enx;
eny_plot=zeros(size(LA));
eny_plot(LA)=eny;
enz_plot=zeros(size(LA));
enz_plot(LA)=enz;

bad=isnan(elx_plot) | isnan(etx_plot) | isnan(enx_plot);
ind_bad=find(bad);
[Ny,Nx,~]=size(LA);
while not(isempty(ind_bad))
    for ii=[1 -1 Ny -Ny Ny*Nx -Ny*Nx]
        upt_bad=not(bad(ind_bad+ii)) & LA(ind_bad+ii);
        upt=ind_bad(upt_bad);

        elx_plot(upt)=elx_plot(upt+ii);
        ely_plot(upt)=ely_plot(upt+ii);
        elz_plot(upt)=elz_plot(upt+ii);

        etx_plot(upt)=etx_plot(upt+ii);
        ety_plot(upt)=ety_plot(upt+ii);
        etz_plot(upt)=etz_plot(upt+ii);

        enx_plot(upt)=enx_plot(upt+ii);
        eny_plot(upt)=eny_plot(upt+ii);
        enz_plot(upt)=enz_plot(upt+ii);

        bad=isnan(elx_plot) | isnan(etx_plot) | isnan(enx_plot);
        ind_bad=find(bad);
    end

end
elx=elx_plot(LA);
ely=ely_plot(LA);
elz=elz_plot(LA);

etx=etx_plot(LA);
ety=ety_plot(LA);
etz=etz_plot(LA);

enx=enx_plot(LA);
eny=eny_plot(LA);
enz=enz_plot(LA);

f.x=elx;
f.y=ely;
f.z=elz;

s.x=etx;
s.y=ety;
s.z=etz;

n.x=enx;
n.y=eny;
n.z=enz;