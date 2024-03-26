function [As] = heartMatrixFib(heartGeom,mask_heart,maskF,fx,fy,fz,D,ratio,res_mm)
%%CREATE ANISOTROPIC DIFFUSION MATRIX FROM INPUTS:%%%%%%%%%%%%%%%%%%%%%%%%%
%heartGeom: phase field voxel matrix of the heart geometry 
%       NOTE: values must be thresholded or you are dividing by 0
%mask_heart: mask of computational domain (where heartGeom > minimal threshold)
%maskF: mask of fibrotic voxels - deprecated-
%fx,fy,fz: x,y,z components of the normalized cardiac fibers
%D: tissue diffusivity 
%ratio: anisotropy ratio
%res_mm: resolution in mm of the input voxel volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_in=mask_heart & ~maskF;
ind = find(ind_in); %computational domain
notmask_heart = not(ind_in); %outside computational domain

%% phase field term [Fenton et al.]
[Ny,Nx,Nz] = size(heartGeom);
nel=Nx*Ny*Nz;
dx = res_mm/10; %resolution in cm
dx2 = dx^2;

%% NOTE:
% Mx,My,Mz assume no flux boundary conditions at the interface of
% mask_heart, resulting in different gradient at the edges when computed considering the
% whole phase field domain. see end of file. No flux should be correct
% because the phase field is assumed to be stationary at the cut threshold

Mx = gradMat(mask_heart,'x',dx);
My = gradMat(mask_heart,'y',dx);
Mz = gradMat(mask_heart,'z',dx);
Mx=Mx(ind_in(mask_heart),:);
My=My(ind_in(mask_heart),:);
Mz=Mz(ind_in(mask_heart),:);

%heartGeom(ind) can not be zero (outside computational domain)
kpx = (Mx*heartGeom(mask_heart))./(heartGeom(ind));
kpy = (My*heartGeom(mask_heart))./(heartGeom(ind));
kpz = (Mz*heartGeom(mask_heart))./(heartGeom(ind));

Dort = D/ratio;
Dpar = D; 
fx = fx(mask_heart);
fy = fy(mask_heart);
fz = fz(mask_heart);
dxfx = Mx*fx;
dyfx = My*fx;
dzfx = Mz*fx;

dxfy = Mx*fy;
dyfy = My*fy;
dzfy = Mz*fy;

dxfz = Mx*fz;
dyfz = My*fz;
dzfz = Mz*fz;
clear Mx My Mz
fx=fx(ind_in(mask_heart));
fy=fy(ind_in(mask_heart));
fz=fz(ind_in(mask_heart));
%% matrix coefficients

%tensor components
d11 = Dort + (Dpar - Dort)*fx.^2;
d12 = (Dpar - Dort)*fx.*fy;
d13 = (Dpar - Dort)*fx.*fz;

d21 = (Dpar - Dort)*fy.*fx;
d22 = Dort + (Dpar - Dort)*fy.^2;
d23 = (Dpar - Dort)*fy.*fz;

d31 = (Dpar - Dort)*fz.*fx;
d32 = (Dpar - Dort)*fz.*fy;
d33 = Dort + (Dpar - Dort)*fz.^2;

%cross derivatives
d12_x1 = (Dpar - Dort)*(fx.*dxfy + fy.*dxfx);
d22_x2 = (Dpar - Dort)*(fy.*dyfy + fy.*dyfy);
d32_x3 = (Dpar - Dort)*(fz.*dzfy + fy.*dzfz);

d11_x1 = (Dpar - Dort)*(fx.*dxfx + fx.*dxfx);
d21_x2 = (Dpar - Dort)*(fy.*dyfx + fx.*dyfy);
d31_x3 = (Dpar - Dort)*(fz.*dzfx + fx.*dzfz);

d13_x1 = (Dpar - Dort)*(fx.*dxfz + fz.*dxfx);
d23_x2 = (Dpar - Dort)*(fy.*dyfz + fz.*dyfy);
d33_x3 = (Dpar - Dort)*(fz.*dzfz + fz.*dzfz);
clear dxfy dyfy dzfy dxfx dyfx dzfx dxfz dyfz dzfz

%% position indexes - 6 near

indi = 1:numel(ind);

%points near the edge of computational box
ind_nodx = logical(mod(ind-1,Nx*Ny)>=Ny*(Nx-1)); % pts with comp. box on the right
ind_nosx = logical(mod(ind-1,Nx*Ny)<Ny); 
ind_nodw = logical(mod(ind,Ny)==0);
ind_noup = logical(mod(ind-1,Ny)==0);
ind_nofr = logical(ind>nel-Nx*Ny);
ind_nobk = logical(ind<= Nx*Ny);
%points near the edge of computational mask
%%NOTE: MUST CHECK INDEXERROR
indout = numel(mask_heart);
ind_nodxout = notmask_heart(min(ind + Ny,indout)); %pts with border on the right
ind_nosxout = notmask_heart(max(ind - Ny,1));
ind_nodwout = notmask_heart(min(ind + 1,indout));
ind_noupout = notmask_heart(max(ind - 1,1));
ind_nofrout = notmask_heart(min(ind + Nx*Ny,indout));
ind_nobkout = notmask_heart(max(ind - Nx*Ny,1));

%merge position indexes
ind_nodx = or(ind_nodx,ind_nodxout);
ind_nosx = or(ind_nosx,ind_nosxout);
ind_nodw = or(ind_nodw,ind_nodwout);
ind_noup = or(ind_noup,ind_noupout);
ind_nofr = or(ind_nofr,ind_nofrout);
ind_nobk = or(ind_nobk,ind_nobkout);
%% position indexes - 12 near
%points near to computational box
ind_nodxdw = notmask_heart(min(ind + Ny+1,indout)); %check mathematical formalism
%points which lack a dx or sx node must be removed as well. 

ind_nodxup = notmask_heart(min(ind + Ny-1,indout));
ind_nodxfr = notmask_heart(min(ind + Ny + Nx*Ny,indout));
ind_nodxbk = notmask_heart(max(ind + Ny - Nx*Ny,1));

ind_nosxdw = notmask_heart(max(ind + 1-Ny,1));
ind_nosxup = notmask_heart(max(ind - 1-Ny,1));
ind_nosxfr = notmask_heart(min(ind - Ny + Nx*Ny,indout));
ind_nosxbk = notmask_heart(max(ind - Ny - Nx*Ny,1));

ind_nobkup = notmask_heart(max(ind - 1 - Nx*Ny,1));
ind_nobkdw = notmask_heart(max(ind + 1 - Nx*Ny,1));
ind_nofrup = notmask_heart(min(ind - 1 + Nx*Ny,indout));
ind_nofrdw = notmask_heart(min(ind + 1 + Nx*Ny,indout));

%% matrix initialization
diagel = zeros(numel(ind),1);
ii = zeros(7*numel(ind),1); 
jj = zeros(7*numel(ind),1);
vvA = zeros(7*numel(ind),1);
%% first/second derivatives: y 
%note: derivative along y axis is along matrix columns (next is +1 row)
% +1 row is +1 in linear index (MATLAB documentation).

ind1=not(ind_nodw); %points on the bottom are excluded
ii(1:nnz(ind1))=indi(ind1); %row position in the matrix, point at (x,y,z)
jj_def = ind(ind1) +1; %column position, point at (x,y + dy,z)
jj(1:nnz(ind1))=jj_def; 
tmp=kpx.*d12/2/dx + kpy.*d22/2/dx + kpz.*d32/2/dx...
    +d12_x1/2/dx + d22_x2/2/dx + d32_x3/2/dx...
    + d22/dx2; 

vvA(1:nnz(ind1)) = tmp((ind1)); %coefficients of the first (plus second) derivative in y
diagel(indi(ind_nodw)) = tmp((ind_nodw)); %points on the bottom are replaced with themselves

%note: derivative along y axis is along matrix columns (previous is -1 row)
ind_sofar = nnz(ii);
ind1=not(ind_noup); %points on top are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);
jj_def = ind(ind1) -1;
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 
tmp =  -kpx.*d12/2/dx - kpy.*d22/2/dx - kpz.*d32/2/dx...
    -d12_x1/2/dx - d22_x2/2/dx - d32_x3/2/dx...
    + d22/dx2;
vvA((ind_sofar+1):(ind_sofar + nnz(ind1)))=tmp((ind1)); 

%points on the top are replaced with themselves
diagel(indi(ind_noup)) = diagel(indi(ind_noup)) + tmp((ind_noup)); 
%coefficient of the second derivative in y
diagel = diagel - 2*d22/dx2;
%% first/second derivatives: x 
%note: derivative along x axis is along matrix rows (next is +1 column)
% +1 column is +Ny in linear index (MATLAB documentation).
ind_sofar = nnz(ii);
ind1=not(ind_nodx); %points at the right border are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1); %row position in the matrix, point at (x,y,z)
jj_def = ind(ind1) + Ny; %column position, point at (x + dx,y,z)
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def;

tmp = kpx.*d11/2/dx + kpy.*d21/2/dx + kpz.*d31/2/dx...
    + d11_x1/2/dx + d21_x2/2/dx + d31_x3/2/dx ...
    + d11/dx2;
vvA((ind_sofar+1):(ind_sofar + nnz(ind1)))=tmp((ind1)); %coefficients of the first (plus second) derivative in x

%points on the border are replaced with themselves
diagel(indi(ind_nodx)) = diagel(indi(ind_nodx)) + tmp((ind_nodx));

%note: derivative along x axis is along matrix rows (previous is -1 column)
ind_sofar = nnz(ii);
ind1=not(ind_nosx); %points at the left border are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);
jj_def = ind(ind1) -Ny;
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 

tmp = -kpx.*d11/2/dx - kpy.*d21/2/dx - kpz.*d31/2/dx...
    - d11_x1/2/dx - d21_x2/2/dx - d31_x3/2/dx...
    + d11/dx2;
vvA((ind_sofar+1):(ind_sofar + nnz(ind1)))=tmp((ind1));

%points on the border are replaced with themselves
diagel(indi(ind_nosx)) = diagel(indi(ind_nosx)) + tmp((ind_nosx));
%coefficient of the second derivative in x
diagel = diagel - 2*d11/dx2;


%% first/second derivatives: z 
%note: derivative along z axis is along matrix 3rd dimension (next is +1 in 3rd dim.)
% +1 in 3rd dimension is +Ny*Nx in linear index (MATLAB documentation).
ind_sofar = nnz(ii);
ind1=not(ind_nofr);%points at the front border are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);%row position in the matrix, point at (x,y,z)
jj_def = ind(ind1) + Ny*Nx; %column position, point at (x,y,z + dz)
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 

tmp = kpx.*d13/2/dx + kpy.*d23/2/dx + kpz.*d33/2/dx...
    + d13_x1/2/dx + d23_x2/2/dx + d33_x3/2/dx ...
    + d33/dx2;
vvA((ind_sofar+1):(ind_sofar + nnz(ind1)))=tmp((ind1)); %coefficients of the first (plus second) derivative in z

%points on the border are replaced with themselves
diagel(indi(ind_nofr)) = diagel(indi(ind_nofr)) + tmp((ind_nofr));

%note: derivative along z axis is along matrix 3rd dimension (prev. is -1 in 3rd dim.)
ind_sofar = nnz(ii);
ind1=not(ind_nobk); %points at the back border are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);
jj_def = ind(ind1) - Nx*Ny;
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 
tmp = -kpx.*d13/2/dx - kpy.*d23/2/dx - kpz.*d33/2/dx...
    - d13_x1/2/dx - d23_x2/2/dx - d33_x3/2/dx...
    + d33/dx2;
vvA((ind_sofar+1):(ind_sofar + nnz(ind1)))=tmp((ind1));

%points on the border are replaced with themselves
diagel(indi(ind_nobk)) = diagel(indi(ind_nobk)) + tmp((ind_nobk));
%coefficient of the second derivative in z
diagel = diagel - 2*d33/dx2;

%add coefficients in the diagonal of the matrix
ind_sofar = nnz(ii);
ii((ind_sofar+1):(ind_sofar + numel(ind))) = indi;
jj((ind_sofar+1):(ind_sofar + numel(ind))) = ind;
vvA((ind_sofar+1):(ind_sofar + numel(ind))) = diagel;

%% second order partial derivatives
%%to compute partial derivatives, points from the 12 neighborhood are
%%needed (e.g. f(x,y+dy,z+dz) )
iic = zeros(12*numel(ind),1); 
jjc = zeros(12*numel(ind),1);
vvAc = zeros(12*numel(ind),1);
diagelc = zeros(numel(indi),1);

%% partial derivatives dxdy 
% (x + dx,y + dy,z)
% points that are excluded ( x+dx or y+dy not defined)
%points which do not have both down and up are excluded also since
%approximation is first order (4 neighbors)
ind1=and( or(not(ind_nodw),not(ind_nodx)), not(ind_nodxdw));

jj_def = ind(ind1)+Ny+1; %column position, point at (x+dx,y+dy,z)
ii_def = indi(ind1); %row position in the matrix, point at (x,y,z)
iic(1:numel(ii_def))= ii_def; 
jjc(1:numel(jj_def))= jj_def; 
tmp = d12/4/dx2...
    +d21/4/dx2;
vvAc(1:numel(ii_def)) = tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = tmp((not(ind1)));

ind_sofar = nnz(iic);
% (x+dx, y-dy, z) 
% points that are excluded
ind1= and( or(not(ind_noup),not(ind_nodx)), not(ind_nodxup)); 
jj_def = ind(ind1)+Ny-1; %column position, point at (x+dx,y-dy,z)
ii_def = indi(ind1); %row position in the matrix, point at (x,y,z)
iic((ind_sofar+1):(ind_sofar + numel(ii_def)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + numel(jj_def)))= jj_def; 
tmp =  -d12/4/dx2...
    -d21/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + numel(ii_def))) = tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

ind_sofar = nnz(iic);
% (x-dx, y+dy, z) 
% points that are excluded
ind1=and(or(not(ind_nodw),not(ind_nosx)), not(ind_nosxdw)); 
jj_def = ind(ind1)-Ny+1; %column position, point at (x-dx,y+dy,z)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = -d12/4/dx2...
    -d21/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

ind_sofar = nnz(iic);
% (x-dx, y-dy, z) 
% points that are excluded
ind1=and(or(not(ind_noup),not(ind_nosx)),not(ind_nosxup)); 
jj_def = ind(ind1)-Ny-1; %column position, point at (x-dx,y-dy,z)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = +d12/4/dx2...
    +d21/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

%% partial derivatives dxdz 
ind_sofar = nnz(iic);
% (x + dx,y,z + dz)
% points that are excluded ( x+dx or z+dz not defined)
ind1 = and(or(not(ind_nodx),not(ind_nofr)),not(ind_nodxfr));
jj_def = ind(ind1) + Ny*Nx + Ny; %column position, point at (x+dx,y,z+dz)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = d31/4/dx2 ...
    + d13/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

ind_sofar = nnz(iic);
% (x + dx,y,z - dz)
% points that are excluded 
ind1 = and(or(not(ind_nodx),not(ind_nobk)),not(ind_nodxbk));
jj_def = ind(ind1) - Ny*Nx + Ny; %column position, point at (x+dx,y,z-dz)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = -d31/4/dx2 ...
    - d13/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

ind_sofar = nnz(iic);
% (x - dx,y,z + dz)
% points that are excluded 
ind1 = and(or(not(ind_nosx),not(ind_nofr)),not(ind_nosxfr));
jj_def = ind(ind1) + Ny*Nx - Ny; %column position, point at (x-dx,y,z+dz)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = -d31/4/dx2...
    - d13/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

ind_sofar = nnz(iic);
% (x - dx,y,z - dz)
% points that are excluded 
ind1 = and(or(not(ind_nosx),not(ind_nobk)),not(ind_nosxbk));
jj_def = ind(ind1) - Ny*Nx - Ny; %column position, point at (x-dx,y,z-dz)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = +d31/4/dx2 ...
    + d13/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));
%% partial derivatives dydz 
ind_sofar = nnz(iic);
% (x,y + dy,z + dz)
% points that are excluded
ind1 = and(or(not(ind_nodw),not(ind_nofr)),not(ind_nofrdw));
jj_def = ind(ind1) + Ny*Nx + 1; %column position, point at (x,y+dy,z+dz)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = +d32/4/dx2 ...
    +d23/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

ind_sofar = nnz(iic);
% (x,y + dy,z - dz)
% points that are excluded
ind1 = and(or(not(ind_nodw),not(ind_nobk)),not(ind_nobkdw));
jj_def = ind(ind1) - Ny*Nx + 1; %column position, point at (x,y+dy,z-dz)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = -d32/4/dx2 ...
    -d23/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

ind_sofar = nnz(iic);
% (x,y - dy,z + dz)
% points that are excluded
ind1 = and(or(not(ind_noup),not(ind_nofr)),not(ind_nofrup));
jj_def = ind(ind1) + Ny*Nx - 1; %column position, point at (x,y-dy,z+dz)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = -d32/4/dx2 ...
    - d23/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

ind_sofar = nnz(iic);
% (x,y - dy,z - dz)
% points that are excluded
ind1 = and(or(not(ind_noup),not(ind_nobk)),not(ind_nobkup));
jj_def = ind(ind1) - Ny*Nx - 1; %column position, point at (x,y-dy,z-dz)
ii_def = indi(ind1);
iic((ind_sofar+1):(ind_sofar + nnz(ind1)))= ii_def; 
jjc((ind_sofar+1):(ind_sofar + nnz(ind1)))= jj_def; 
tmp = +d32/4/dx2 ...
    + d23/4/dx2;
vvAc((ind_sofar+1):(ind_sofar + nnz(ind1))) =  tmp((ind1));
%previously excluded points are replaced with themselves
diagelc(indi(not(ind1))) = diagelc(indi(not(ind1)))+tmp((not(ind1)));

%add coefficients on the diagonal of the matrix 
ind_sofar = nnz(iic);
iic((ind_sofar+1):(ind_sofar + numel(ind))) = indi;
jjc((ind_sofar+1):(ind_sofar + numel(ind))) = ind;
vvAc((ind_sofar+1):(ind_sofar + numel(ind))) = diagelc;
%% matrix construction 

%%cut unused vectror space
ind_sofar = nnz(ii);
ii = ii(1:ind_sofar);
jj = jj(1:ind_sofar);
vvA = vvA(1:ind_sofar);
ind_sofar = nnz(iic);
iic = iic(1:ind_sofar);
jjc = jjc(1:ind_sofar);
vvAc = vvAc(1:ind_sofar);
%%allocation
As = sparse(ii, jj, vvA, numel(ind), nel);
Asc = sparse(iic,jjc,vvAc,numel(ind),nel);
As = (As + Asc);
clear Asc
%%a lot of outer points are included in the columns of As.
%%since all outer points were excluded from computation,
%%we only keep rows corresponding to points inside the domain (pointed by ind)
%%the columns pointed by ind are ordered as the rows of As, since they are
%%extracted in the same way. 
clear fx fy fz d11 d12 d13 d21 d22 d23 d31 d32 d33 dxfx dxfy dxfz dyfx dyfy dyfz ...
    dzfx dzfy dzfz d12_x1 d22_x2 d32_x3 d11_x1 d21_x2 d31_x3 d13_x1 d23_x2 d33_x3 
clear diagel diagelc ii iic jj jjc vvA vvAc ind_nodx ind_nosx ind_nodw ...
    ind_noup ind_nofr ind_nobk ind_nodxF ind_nosxF ind_nodwF ind_noupF ...
    ind_nofrF ind_nobkF ind_nodxout ind_nosxout ind_nodwout ind_noupout...
    ind_nofrout ind_nobkout




%%remove unused columns 
As = As(indi,ind);
end
