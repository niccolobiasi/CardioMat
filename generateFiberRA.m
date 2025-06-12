function [f,s,n]=generateFiberRA(RA, fields_RA, par)
% generateFiberLA(LA,fields_LA,par) generate fiber orientation for the
% voxelized right-atrial geometry RA by using a rule-based method (see
% Piersanti et al. 2021). The function uses the distance fields fields_RA
% as computed by computeFieldsRA. Optionally, a par structure can be input
% defining the parameter values for bundles selection with subfields:
%   - tau_tv (default:0.95)
%   - tau_icv (default: 0.95)
%   - tau_scv (default: 0.05)
%   - tau_ct_p (default: -0.05)
%   - tau_ct_m (default: -0.13)
%   - tau_ib (default: -0.1)
%   - tau_ras (default: 0.13)
%   - tau_raw (default: 0.55)

if isfield(par,'tau_tv')
    tau_tv=par.tau_tv;
else
    tau_tv=0.95;
end

if isfield(par,'tau_icv')
    tau_icv=par.tau_icv;
else
    tau_icv=0.95;
end

if isfield(par,'tau_scv')
    tau_scv=par.tau_scv;
else
    tau_scv=0.05;
end

if isfield(par,'tau_ct_p')
    tau_ct_p=par.tau_ct_p;
else
    tau_ct_p=-0.05;
end

if isfield(par,'tau_ct_m')
    tau_ct_m=par.tau_ct_m;
else
    tau_ct_m=-0.13;
end

if isfield(par,'tau_ib')
    tau_ib=par.tau_ib;
else
    tau_ib=0.1;
end

if isfield(par,'tau_ras')
    tau_ras=par.tau_ras;
else
    tau_ras=0.13;
end

if isfield(par,'tau_raw')
    tau_raw=par.tau_raw;
else
    tau_raw=0.55;
end

phi_r=fields_RA{1};
psi_ab_ra=fields_RA{2};
psi_v_ra=fields_RA{3};
psi_r_ra=fields_RA{4};
psi_w_ra=fields_RA{5};




S=gradMat(RA,'x');
tmx=S*phi_r(RA);
r_r_x=S*psi_r_ra(RA);
v_r_x=S*psi_v_ra(RA);
ab_r_x=S*psi_ab_ra(RA);
w_r_x=S*psi_w_ra(RA);

S=gradMat(RA,'y');
tmy=S*phi_r(RA);
r_r_y=S*psi_r_ra(RA);
v_r_y=S*psi_v_ra(RA);
ab_r_y=S*psi_ab_ra(RA);
w_r_y=S*psi_w_ra(RA);

S=gradMat(RA,'z');
tmz=S*phi_r(RA);
r_r_z=S*psi_r_ra(RA);
v_r_z=S*psi_v_ra(RA);
ab_r_z=S*psi_ab_ra(RA);
w_r_z=S*psi_w_ra(RA);

norm_et=sqrt(tmx.^2+tmy.^2+tmz.^2);
etx=tmx./(norm_et);
ety=tmy./(norm_et);
etz=tmz./(norm_et);

kx=zeros(size(etx));
ky=zeros(size(etx));
kz=zeros(size(etx));

ind1=psi_r_ra(RA)>=tau_tv;

kx(ind1)=r_r_x(ind1);
ky(ind1)=r_r_y(ind1);
kz(ind1)=r_r_z(ind1);

ind2=psi_r_ra(RA)<tau_raw;

ind3=psi_w_ra(RA)>=tau_ct_m & psi_w_ra(RA)<=tau_ct_p;

cind=~ind1 & ind2 & ind3;

kx(cind)=w_r_x(cind);
ky(cind)=w_r_y(cind);
kz(cind)=w_r_z(cind);

ind4=psi_w_ra(RA)<tau_ct_m;

ind5=(psi_v_ra(RA)>=tau_icv | psi_v_ra(RA)<=tau_scv);

cind=~ind1 & ind2 & ind4 & ind5;
kx(cind)=v_r_x(cind);
ky(cind)=v_r_y(cind);
kz(cind)=v_r_z(cind);

cind=~ind1 & ind2 & ind4 & ~ind5;
kx(cind)=ab_r_x(cind);
ky(cind)=ab_r_y(cind);
kz(cind)=ab_r_z(cind);


cind=~ind1 & ind2 & ~ind3 & ~ind4 & ind5;
kx(cind)=v_r_x(cind);
ky(cind)=v_r_y(cind);
kz(cind)=v_r_z(cind);

ind6= psi_w_ra(RA)<=tau_ib;


cind=~ind1 & ind2 & ~ind3 & ~ind4 & ~ind5 & ind6;
kx(cind)=v_r_x(cind);
ky(cind)=v_r_y(cind);
kz(cind)=v_r_z(cind);


ind7=psi_w_ra(RA)>=tau_ras;

cind=~ind1 & ind2 & ~ind3 & ~ind4 & ~ind5 & ~ind6 & ind7;
kx(cind)=r_r_x(cind);
ky(cind)=r_r_y(cind);
kz(cind)=r_r_z(cind);

cind=~ind1 & ind2 & ~ind3 & ~ind4 & ~ind5 & ~ind6 & ~ind7;
kx(cind)=w_r_x(cind);
ky(cind)=w_r_y(cind);
kz(cind)=w_r_z(cind);


cind=~ind1 & ~ind2 & ind5;
kx(cind)=v_r_x(cind);
ky(cind)=v_r_y(cind);
kz(cind)=v_r_z(cind);

ind8=psi_w_ra(RA)>=0;

cind=~ind1 & ~ind2 & ~ind5 & ind8;
kx(cind)=r_r_x(cind);
ky(cind)=r_r_y(cind);
kz(cind)=r_r_z(cind);

cind=~ind1 & ~ind2 & ~ind5 & ~ind8;
kx(cind)=ab_r_x(cind);
ky(cind)=ab_r_y(cind);
kz(cind)=ab_r_z(cind);

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

elx_plot=zeros(size(RA));
elx_plot(RA)=elx;
ely_plot=zeros(size(RA));
ely_plot(RA)=ely;
elz_plot=zeros(size(RA));
elz_plot(RA)=elz;

etx_plot=zeros(size(RA));
etx_plot(RA)=etx;
ety_plot=zeros(size(RA));
ety_plot(RA)=ety;
etz_plot=zeros(size(RA));
etz_plot(RA)=etz;

enx_plot=zeros(size(RA));
enx_plot(RA)=enx;
eny_plot=zeros(size(RA));
eny_plot(RA)=eny;
enz_plot=zeros(size(RA));
enz_plot(RA)=enz;

bad=isnan(elx_plot) | isnan(etx_plot) | isnan(enx_plot);
ind_bad=find(bad);
[Ny,Nx,~]=size(RA);
while not(isempty(ind_bad))
    for ii=[1 -1 Ny -Ny Ny*Nx -Ny*Nx]
        upt_bad=not(bad(ind_bad+ii)) & RA(ind_bad+ii);
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
elx=elx_plot(RA);
ely=ely_plot(RA);
elz=elz_plot(RA);

etx=etx_plot(RA);
ety=ety_plot(RA);
etz=etz_plot(RA);

enx=enx_plot(RA);
eny=eny_plot(RA);
enz=enz_plot(RA);

f.x=elx;
f.y=ely;
f.z=elz;

s.x=etx;
s.y=ety;
s.z=etz;

n.x=enx;
n.y=eny;
n.z=enz;