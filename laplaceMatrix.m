function S=laplaceMatrix(Mat)
% S=laplaceMatrix(siz) returns the finite difference matrix for the
% computation of laplacian operator of a generic field, such that 
% S*field/dx^2 is the laplacian of field.
[Ny,Nx,Nz]=size(Mat);
nel=Nx*Ny*Nz;

ind=find(Mat);

indi = 1:numel(ind);

%points near the edge of computational box
ind_nodx = logical(mod(ind-1,Nx*Ny)>=Ny*(Nx-1)); % pts with comp. box on the right
ind_nosx = logical(mod(ind-1,Nx*Ny)<Ny); 
ind_nodw = logical(mod(ind,Ny)==0);
ind_noup = logical(mod(ind-1,Ny)==0);
ind_nofr = logical(ind>nel-Nx*Ny);
ind_nobk = logical(ind<= Nx*Ny);
%points near the edge of computational mask
notmask_heart=~Mat;
ind_nodxout = notmask_heart(min(ind + Ny,nel)); %pts with border on the right
ind_nosxout = notmask_heart(max(ind - Ny,1));
ind_nodwout = notmask_heart(min(ind + 1,nel));
ind_noupout = notmask_heart(max(ind - 1,1));
ind_nofrout = notmask_heart(min(ind + Nx*Ny,nel));
ind_nobkout = notmask_heart(max(ind - Nx*Ny,1));

%merge position indexes
ind_nodx = or(ind_nodx,ind_nodxout);
ind_nosx = or(ind_nosx,ind_nosxout);
ind_nodw = or(ind_nodw,ind_nodwout);
ind_noup = or(ind_noup,ind_noupout);
ind_nofr = or(ind_nofr,ind_nofrout);
ind_nobk = or(ind_nobk,ind_nobkout);

ii = zeros(7*numel(ind),1); 
jj = zeros(7*numel(ind),1);

ind1=not(ind_nodw); %points on the bottom are excluded
ii(1:nnz(ind1))=indi(ind1); %row position in the matrix, point at (x,y,z)
jj_def = ind(ind1) +1; %column position, point at (x,y + dy,z)
jj(1:nnz(ind1))=jj_def; 

ind_sofar = nnz(ii);
ind1=not(ind_noup); %points on top are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);
jj_def = ind(ind1) -1;
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def;

ind_sofar = nnz(ii);
ind1=not(ind_nodx); %points at the right border are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1); %row position in the matrix, point at (x,y,z)
jj_def = ind(ind1) + Ny; %column position, point at (x + dx,y,z)
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def;

ind_sofar = nnz(ii);
ind1=not(ind_nosx); %points at the left border are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);
jj_def = ind(ind1) -Ny;
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 


ind_sofar = nnz(ii);
ind1=not(ind_nofr);%points at the front border are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);%row position in the matrix, point at (x,y,z)
jj_def = ind(ind1) + Ny*Nx; %column position, point at (x,y,z + dz)
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 

ind_sofar = nnz(ii);
ind1=not(ind_nobk); %points at the back border are excluded
ii((ind_sofar+1):(ind_sofar + nnz(ind1)))= indi(ind1);
jj_def = ind(ind1) - Nx*Ny;
jj((ind_sofar+1):(ind_sofar + nnz(ind1)))=jj_def; 

ind_sofar = nnz(ii);
ii = ii(1:ind_sofar);
jj = jj(1:ind_sofar);

S=sparse(ii,jj,ones(length(ii),1),numel(ind),nel,length(ii)+nel);
S = S(indi,ind);
S=S-diag(sum(S,2));
end