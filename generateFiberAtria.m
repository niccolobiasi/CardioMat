function [f,s,n]=generateFiberAtria(VoxelMat,LA,RA,fields_LA, fields_RA,par_la,par_ra)
% generateFiberAtria(VoxelMat,LA,RA,fields_LA, fields_RA,par_la,par_ra) 
% generate fiber orientation for the voxelized biatrial geometry VoxelMat
% by using a rule-based method (see Piersanti et al. 2021).
% The function uses the atria separation (LA, and RA), and the distance 
% fields (fields_LA and fields_RA) as computed by computeDistanceFields. 
% Optionally, for each atrium a par structure can be input defining 
% the parameter values for bundles selection. For further details see
% generateFiberLA and generateFiberRA documentation.

if nargin<7
    par_ra=[];
end
if nargin<6
    par_la=[];
end
[f_la,s_la,n_la]=generateFiberLA(LA,fields_LA,par_la);

[f_ra,s_ra,n_ra]=generateFiberLA(RA,fields_RA,par_ra);


f.x=zeros(nnz(VoxelMat),1);
f.y=zeros(nnz(VoxelMat),1);
f.z=zeros(nnz(VoxelMat),1);

f.x(LA(VoxelMat))=f_la.x;
f.y(LA(VoxelMat))=f_la.y;
f.z(LA(VoxelMat))=f_la.z;

f.x(RA(VoxelMat))=f_ra.x;
f.y(RA(VoxelMat))=f_ra.y;
f.z(RA(VoxelMat))=f_ra.z;

n.x=zeros(nnz(VoxelMat),1);
n.y=zeros(nnz(VoxelMat),1);
n.z=zeros(nnz(VoxelMat),1);

n.x(LA(VoxelMat))=n_la.x;
n.y(LA(VoxelMat))=n_la.y;
n.z(LA(VoxelMat))=n_la.z;

n.x(RA(VoxelMat))=n_ra.x;
n.y(RA(VoxelMat))=n_ra.y;
n.z(RA(VoxelMat))=n_ra.z;

s.x=zeros(nnz(VoxelMat),1);
s.y=zeros(nnz(VoxelMat),1);
s.z=zeros(nnz(VoxelMat),1);

s.x(LA(VoxelMat))=s_la.x;
s.y(LA(VoxelMat))=s_la.y;
s.z(LA(VoxelMat))=s_la.z;

s.x(RA(VoxelMat))=s_ra.x;
s.y(RA(VoxelMat))=s_ra.y;
s.z(RA(VoxelMat))=s_ra.z;