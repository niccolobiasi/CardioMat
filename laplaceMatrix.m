function S=laplaceMatrix(siz)
% S=laplaceMatrix(siz) returns the finite difference matrix for the
% computation of laplacian operator of a generic field, such that 
% S*field/dx^2 is the laplacian of field.
Ny=siz(1);
Nx=siz(2);
Nz=siz(3);
nel=Nx*Ny*Nz;

ind=(1:nel)';

forw_dy=not(mod(ind,Ny)==0);
back_dy=not(mod(ind-1,Ny)==0) & ind>1;
forw_dx=mod(ind-1,Nx*Ny)<Ny*(Nx-1);
back_dx=mod(ind-1,Nx*Ny)>=Ny;
forw_dz=ind<=nel-Nx*Ny;
back_dz=ind>Nx*Ny;

ind1=ind(forw_dy);
ii=ind1;
jj=ind1+1;

ind1=ind(back_dy);
ii=[ii; ind1];
jj=[jj; ind1-1];

ind1=ind(forw_dx);
ii=[ii; ind1];
jj=[jj; ind1+Ny];

ind1=ind(back_dx);
ii=[ii; ind1];
jj=[jj; ind1-Ny];

ind1=ind(forw_dz);
ii=[ii; ind1];
jj=[jj; ind1+Nx*Ny];

ind1=ind(back_dz);
ii=[ii; ind1];
jj=[jj; ind1-Nx*Ny];

S=sparse(ii,jj,ones(length(ii),1),nel,nel,length(ii)+nel);
S=S-diag(sum(S,2));
end