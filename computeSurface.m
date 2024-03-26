function [FV,extInd]=computeSurface(VoxelMat,res,minBound)
% computeSurface finds the surface voxels in the geometry as inside voxels
% with adiacent voxels outside the geometry.
%[FV,extInd]=computeSurface(VoxelMat,res) returns a structure FV containing
%two subfields (vertices and faces) and a 1D array extInd. extInd contains
%the voxel indices corresponding to each face.
% Note that this function assumes that true values of VoxelMat are not on
% the bounadries of the computational domain
%
if nargin<3
    Xmin=0;
    Ymin=0;
    Zmin=0;
else
    Xmin=minBound(1);
    Ymin=minBound(2);
    Zmin=minBound(3);
end
dx2=res/2;
[Ny,Nx,~]=size(VoxelMat);
ind_in=find(VoxelMat);
noup=~VoxelMat(ind_in+1);
nodw=~VoxelMat(ind_in-1);
nodx=~VoxelMat(ind_in+Ny);
nosx=~VoxelMat(ind_in-Ny);
nofw=~VoxelMat(ind_in+Ny*Nx);
nobw=~VoxelMat(ind_in-Ny*Nx);
ind_surf=find(noup | nodw | nodx | nosx | nofw | nobw);

[jy,jx,jz]=ind2sub(size(VoxelMat),ind_in(ind_surf));
jx=(repmat(jx,[1, 8]))';
jx=jx(:);
jy=(repmat(jy,[1, 8]))';
jy=jy(:);
jz=(repmat(jz,[1, 8]))';
jz=jz(:);
FV.vertices=res*[jx jy jz]-dx2+repmat([dx2 dx2 dx2; dx2 -dx2 dx2;...
        dx2 -dx2 -dx2; dx2 dx2 -dx2; -dx2 dx2 dx2; -dx2 -dx2 dx2;...
        -dx2 -dx2 -dx2; -dx2 dx2 -dx2],[length(ind_surf) 1]);

nodx=find(nodx(ind_surf));
FV.faces=repmat([1 2 3 4],[length(nodx) 1])+(nodx-1)*8;
extInd=ind_in(ind_surf(nodx));

nosx=find(nosx(ind_surf));
FV.faces=[FV.faces; repmat([6 5 8 7],[length(nosx) 1])+(nosx-1)*8];
extInd=[extInd; ind_in(ind_surf(nosx))];

noup=find(noup(ind_surf));
FV.faces=[FV.faces; repmat([5 1 4 8],[length(noup) 1])+(noup-1)*8];
extInd=[extInd; ind_in(ind_surf(noup))];

nodw=find(nodw(ind_surf));
FV.faces=[FV.faces; repmat([2 6 7 3],[length(nodw) 1])+(nodw-1)*8];
extInd=[extInd; ind_in(ind_surf(nodw))];

nofw=find(nofw(ind_surf));
FV.faces=[FV.faces; repmat([1 2 6 5],[length(nofw) 1])+(nofw-1)*8];
extInd=[extInd; ind_in(ind_surf(nofw))];

nobw=find(nobw(ind_surf));
FV.faces=[FV.faces; repmat([4 3 7 8],[length(nobw) 1])+(nobw-1)*8];
extInd=[extInd; ind_in(ind_surf(nobw))];

FV.vertices(:,1)=FV.vertices(:,1)+Ymin;
FV.vertices(:,2)=FV.vertices(:,2)+Xmin;
FV.vertices(:,3)=FV.vertices(:,3)+Zmin;

end