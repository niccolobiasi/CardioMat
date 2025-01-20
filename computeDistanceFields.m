function [transatrial, LA, RA, fields_LA, fields_RA, phi]=computeDistanceFields(VoxelMat,tag,res)
% [transatrial, LA, RA, fields_LA, fields_RA, phi]=computeDistanceFields(VoxelMat,tag)
% returns the transatrial coordinate and separates left (LA) and right (RA)
% atria. The function also computes distance fields useful for the 
% assignment of fiber orientation both in the left (fields_LA) and right
% (fields_RA) atria. The biatrial transmural coordinate phi can be output
% as well. The function requires the biatrial geometry VoxelMat and the
% surface tags tag as produced by tagAtria function. The resolution res is
% optional (default: 0.25mm)

if nargin<3
    res=0.25;
end

eik=solveEikonal(VoxelMat,find(tag==8));
eik1=solveEikonal(VoxelMat,find(tag==7));
transatrial=eik./(eik+eik1);
if nargout>1
    RA=reshape(transatrial>0.5,size(VoxelMat)) & VoxelMat;
    LA=reshape(transatrial<=0.5,size(VoxelMat)) & VoxelMat;
end

if nargout>3
    fields_LA=computeFieldsLA(LA,tag,res);
end

if nargout>4
    fields_RA=computeFieldsRA(RA,tag,res);
end

if nargout>5
    phi=0.5*ones(size(VoxelMat));
    phi(LA)=fields_LA{1}(LA);
    phi(RA)=fields_RA{1}(RA);
    phi=phi(:);
end