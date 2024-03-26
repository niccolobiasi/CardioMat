function As=gradMat(VoxelMat,dir,dx)
% S=gradMat(VoxelMat,dir,dx) returns a matrix S such that S*sf is the 
% approximated first derivative of the generic scalar field sf in the 
% direction specified by dir. 
% VOXELMAT:  logical mask where the scalar field is defined.
% DIR: direction along which the gradient is computed (char 'x', 'y', 'z')
% DX: resolution of the input mask


if nargin<3
    dx=1;
end

notmask_heart = not(VoxelMat); %outside computational domain
ind = find(VoxelMat); %computational domain
[Ny,Nx,Nz]=size(VoxelMat);
nel=Nx*Ny*Nz;
%% position indexes - 6 near
indi = 1:numel(ind);

%points near to computational box
ind_nodx = logical(mod(ind-1,Nx*Ny)>=Ny*(Nx-1)); % pts with comp. box on the right
ind_nosx = logical(mod(ind-1,Nx*Ny)<Ny); 
ind_nodw = logical(mod(ind,Ny)==0);
ind_noup = logical(mod(ind-1,Ny)==0);
ind_nofr = logical(ind>nel-Nx*Ny);
ind_nobk = logical(ind<= Nx*Ny);
%points near the mask computational border
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
diagel = zeros(numel(ind),1);
ii = zeros(3*numel(ind),1); 
jj = zeros(3*numel(ind),1);
vvA = zeros(3*numel(ind),1);

if dir =='y'
%% first/second derivatives: y 
%note: derivative along y axis is along matrix columns (next is +1 row)
% +1 row is +1 in linear index (MATLAB documentation).

ind1=not(ind_nodw); %points excluded
ii(1:nnz(ind1))=indi(ind1); %row position in the gradient matrix, points to (x,y,z)
jj_def = ind(ind1) +1; %column position, points to (x,y + dy,z)
jj(1:nnz(ind1))=jj_def; 

tmp=1/2/dx; %coefficient 
vvA(1:nnz(ind1)) = tmp; 
diagel(indi(ind_nodw)) = tmp; %points excluded are replaced with themselves

%note: derivative along y axis is along matrix columns (previous is -1 row)
ind_sofar = nnz(ii);
ind1=not(ind_noup); %points excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);
jj_def = ind(ind1) -1;
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 

tmp =  -1/2/dx;
vvA((ind_sofar+1):(ind_sofar + nnz(ind1)))=tmp; 
diagel(indi(ind_noup)) = diagel(indi(ind_noup)) + tmp; 
elseif dir =='x'
%% first/second derivatives: x 
%note: derivative along x axis is along matrix rows (next is +1 column)
% +1 column is +Ny in linear index (MATLAB documentation).
ind1=not(ind_nodx); %points  excluded
ii(1:nnz(ind1))=indi(ind1); %row position in the gradient matrix, points to (x,y,z)
jj_def = ind(ind1) + Ny; %column position, points to (x + dx,y,z)
jj(1:nnz(ind1))=jj_def; 

tmp = 1/2/dx;
vvA(1:nnz(ind1)) = tmp; 
diagel(indi(ind_nodx)) = tmp; %points excluded are replaced with themselves

%note: derivative along x axis is along matrix rows (previous is -1 column)
ind_sofar = nnz(ii);
ind1=not(ind_nosx); 
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);
jj_def = ind(ind1) -Ny;
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 

tmp = -1/2/dx;
vvA((ind_sofar+1):(ind_sofar + nnz(ind1)))=tmp;
diagel(indi(ind_nosx)) = diagel(indi(ind_nosx)) + tmp;
elseif dir == 'z'
%% first/second derivatives: z 
%note: derivative along z axis is along matrix 3rd dimension (next is +1 in 3rd dim.)
% +1 in 3rd dimension is +Ny*Nx in linear index (MATLAB documentation).
ind1=not(ind_nofr);%points at the front border are excluded
ii(1:nnz(ind1))=indi(ind1); %row position in the gradient matrix, points to (x,y,z)
jj_def = ind(ind1) + Ny*Nx; %column position, point to (x,y,z + dz)
jj(1:nnz(ind1))=jj_def; 

tmp = 1/2/dx;
vvA(1:nnz(ind1)) = tmp;  
diagel(indi(ind_nofr)) = tmp; %points excluded are replaced with themselves

%note: derivative along z axis is along matrix 3rd dimension (prev. is -1 in 3rd dim.)
ind_sofar = nnz(ii);
ind1=not(ind_nobk); 
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);
jj_def = ind(ind1) - Nx*Ny;
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 

tmp = -1/2/dx;
vvA((ind_sofar+1):(ind_sofar + nnz(ind1)))=tmp;
diagel(indi(ind_nobk)) = diagel(indi(ind_nobk)) + tmp;
end

%allocate coefficients on the diagonal
ind_sofar = nnz(ii);
ii((ind_sofar+1):(ind_sofar + numel(ind))) = indi;
jj((ind_sofar+1):(ind_sofar + numel(ind))) = ind;
vvA((ind_sofar+1):(ind_sofar + numel(ind))) = diagel;

%cut unused vector space
ind_sofar = nnz(ii);
ii = ii(1:ind_sofar);
jj = jj(1:ind_sofar);
vvA = vvA(1:ind_sofar);

%%allocation
As = sparse(ii, jj, vvA, numel(ind), nel);
%%a lot of outer points are included in the columns of As.
%%since all outer points were excluded from computation,
%%we only keep rows corresponding to points inside the domain (pointed by ind)
%%the columns pointed by ind are ordered as the rows of As, since they are
%%extracted in the same way. 

%%remove unused columns 
As = As(indi,ind);

end