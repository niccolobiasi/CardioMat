function [LF] = LFvector(heartGeom,mask_heart,maskF,fx,fy,fz,D,ratio,res_mm,pos)
% This function creates a matrix for computation of extracellular potential
% in the position identified by the matrix pos. This code assumes the tissu
% e to be immersed in an infinite homogeneous volume conductor.
% INPUTS
% heartGeom: phase field voxel matrix of the heart geometry
% mask_heart: computational domain (where phi > minimal threshold)
% maskF: mask of fibrotic voxels
% fx,fy,fz: x,y,z components of the normalized cardiac fibers
% D: tissue diffusivity 
% ratio: anisotropy ratio
% res_mm: resolution in mm of the input voxel volume
% pos: Nx3 matrix with each row corresponding to a point where to compute
% the extracellular potential.
% OUTPUT
% LF: a (N x nin) matrix such that 1/(4*pi*sigma_bath)*LF*Vm returns a
% vector with the extracellular potential at the N positions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_in=mask_heart & ~maskF;
ind = find(ind_in);
notmask_heart = not(ind_in); 

%% phase field term [Fenton et al.]
[Ny,Nx,Nz] = size(heartGeom);
nel=Nx*Ny*Nz;
fx = fx(ind); fy = fy(ind); fz = fz(ind);
heartGeom = heartGeom(ind);
%% Diffusivity tensor
Dort =  D/ratio;
Dpar = D; 
d11 = Dort + (Dpar - Dort)*fx.^2;
d12 = (Dpar - Dort)*fx.*fy;
d13 = (Dpar - Dort)*fx.*fz;

d21 = (Dpar - Dort)*fy.*fx;
d22 = Dort + (Dpar - Dort)*fy.^2;
d23 = (Dpar - Dort)*fy.*fz;

d31 = (Dpar - Dort)*fz.*fx;
d32 = (Dpar - Dort)*fz.*fy;
d33 = Dort + (Dpar - Dort)*fz.^2;




%% position indexes - 6 near

%points near to computational box
ind_nodx = logical(mod(ind-1,Nx*Ny)>=Ny*(Nx-1)); % pts with comp. box on the right
ind_nosx = logical(mod(ind-1,Nx*Ny)<Ny); 
ind_nodw = logical(mod(ind,Ny)==0);
ind_noup = logical(mod(ind-1,Ny)==0);
ind_nofr = logical(ind>nel-Nx*Ny);
ind_nobk = logical(ind<= Nx*Ny);
%points near the computational border
ind_nodxout = notmask_heart(ind + Ny); %pts with border on the right
ind_nosxout = notmask_heart(ind - Ny);
ind_nodwout = notmask_heart(ind + 1);
ind_noupout = notmask_heart(ind - 1);
ind_nofrout = notmask_heart(ind + Nx*Ny);
ind_nobkout = notmask_heart(ind - Nx*Ny);

%merge position indexes
ind_nodx = or(ind_nodx,ind_nodxout);
ind_nosx = or(ind_nosx,ind_nosxout);
ind_nodw = or(ind_nodw,ind_nodwout);
ind_noup = or(ind_noup,ind_noupout);
ind_nofr = or(ind_nofr,ind_nofrout);
ind_nobk = or(ind_nobk,ind_nobkout);

%% matrix initialization
[iiq,jjq,kkq] = ind2sub(size(mask_heart),ind);
LF = zeros(size(pos,1),length(ind));
for posit = 1:size(pos,1)

w = zeros(nel,1);
rx = jjq - pos(posit,1)/res_mm; %coord x
ry = iiq - pos(posit,2)/res_mm; %coord y
rz = kkq - pos(posit,3)/res_mm; %coord z
R = sqrt(ry.^2 + rx.^2 + rz.^2).^3;
LFx = (rx./(R+eps)); LFy = (ry./(R+eps)); LFz = (rz./(R+eps));
%% first/second derivatives: y 
%note: derivative along y axis is along matrix columns (next is +1 row)
% +1 row is +1 in linear index (MATLAB documentation).

ind1=not(ind_nodw); %position P
jj_def = ind(ind1) +1; %position, at (x,y + dy,z)
coefs = (LFx.*d12 + LFy.*d22 + LFz.*d32).*heartGeom;
%points with a next
w(jj_def) =  coefs(ind1);
w(ind(ind_nodw)) = w(ind(ind_nodw))+ coefs(ind_nodw);

%note: derivative along y axis is along matrix columns (previous is -1 row)
ind1=not(ind_noup); %points on top are excluded
jj_def = ind(ind1) -1;
%points with a next
w(jj_def) = w(jj_def) - coefs(ind1);
w(ind(ind_noup)) = w(ind(ind_noup)) - coefs(ind_noup);

%% first/second derivatives: x 
%note: derivative along x axis is along matrix rows (next is +1 column)
% +1 column is +Ny in linear index (MATLAB documentation).

ind1=not(ind_nodx); %points at the right border are excluded
jj_def = ind(ind1) + Ny; %column position, point at (x + dx,y,z)
coefs = (LFx.*d11 + LFy.*d21 + LFz.*d31).*heartGeom;

w(jj_def) = w(jj_def) + coefs(ind1);
w(ind(ind_nodx)) = w(ind(ind_nodx)) + coefs(ind_nodx);

ind1=not(ind_nosx); %points at the left border are excluded
jj_def = ind(ind1) -Ny;
w(jj_def) = w(jj_def) - coefs(ind1);
w(ind(ind_nosx)) = w(ind(ind_nosx)) - coefs(ind_nosx);

%% first/second derivatives: z 
%note: derivative along z axis is along matrix 3rd dimension (next is +1 in 3rd dim.)
% +1 in 3rd dimension is +Ny*Nx in linear index (MATLAB documentation).

ind1=not(ind_nofr);%points at the front border are excluded
jj_def = ind(ind1) + Ny*Nx; %column position, point at (x,y,z + dz)
coefs = (LFx.*d13 + LFy.*d23 + LFz.*d33).*heartGeom;
w(jj_def) = w(jj_def) + coefs(ind1);
w(ind(ind_nofr)) = w(ind(ind_nofr)) + coefs(ind_nofr);
%note: derivative along z axis is along matrix 3rd dimension (prev. is -1 in 3rd dim.)
ind1=not(ind_nobk); %points at the back border are excluded
jj_def = ind(ind1) - Nx*Ny;
w(jj_def) = w(jj_def) - coefs(ind1);
w(ind(ind_nobk)) = w(ind(ind_nobk)) - coefs(ind_nobk);


%keep only values in the computational domain
w = w(ind)/2;

LF(posit,:) = w';
end



end

