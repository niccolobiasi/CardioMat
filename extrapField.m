function f_new=extrapField(f,VoxelMat,extInd,mask_phi,max_it)
%f_new=extrapField(f,VoxelMat,extInd) returns the extrapolation in the
%whole computational domain (as defined by the size of the voxelized 
% geometry VoxelMat)of the scalar fields contained in f.
%f is a ninXnFields matrix where nin=nnz(VoxelMat) and nFields is the
%number of fields to extrapolate.
%extInd is a vector containing the indices of the voxels in the surface of
%the voxelized geometry.
%
%f_new=extrapField(f,VoxelMat,extInd, mask_phi) will stop the extrapolation
%when all the voxels marked by the mask_phi matrix have a defined value for
%the fields. The remaining voxels at the stop of the function are left NaN
%values.
%
%f_new=extrapField(f,VoxelMat,extInd, mask_phi,max_it) will stop the
%extrapolation when all the voxels marked by the mask_phi matrix have a 
% defined value for the fields or the maximum number of iteration is
% achieved. Default maximum number of iteration is infinite.

if size(f,1)~=nnz(VoxelMat)
    error('The length of vector field must be equal to the number of internal voxels')
end

nel=numel(VoxelMat);
if nargin<4
    mask_phi=true(nel,1);
else
    if isempty(mask_phi)
        mask_phi=true(nel,1);
    else
        mask_phi=mask_phi(:);
    end
end

if nargin<5
    max_it=Inf;
end

f_new=NaN([nel,size(f,2)]);
f_new(VoxelMat(:),:)=f;
[Ny,Nx,~]=size(VoxelMat);
seed=unique(extInd);
it=0;
while nnz(isnan(f_new(mask_phi,1)))>0 && it<max_it
    if isempty(seed)
        disp('there are unreacheable voxels')
        break;
    end

    new_seed=[];
    
    sel_seed=seed(mod(seed,Ny)>0);
    fw_y=sel_seed+1;
    ToUpdate=isnan(f_new(fw_y,1));
    new_seed=[new_seed; fw_y(ToUpdate)];
    f_new(fw_y(ToUpdate),:)=f_new(sel_seed(ToUpdate),:);

    sel_seed=seed(mod(seed-1,Ny)>0 & seed>1);
    bw_y=sel_seed-1;
    ToUpdate=isnan(f_new(bw_y,1));
    new_seed=[new_seed; bw_y(ToUpdate)];
    f_new(bw_y(ToUpdate),:)=f_new(sel_seed(ToUpdate),:);

    sel_seed=seed(mod(seed-1,Nx*Ny)<Ny*(Nx-1));
    fw_x=sel_seed+Ny;
    ToUpdate=isnan(f_new(fw_x,1));
    new_seed=[new_seed; fw_x(ToUpdate)];
    f_new(fw_x(ToUpdate),:)=f_new(sel_seed(ToUpdate),:);

    sel_seed=seed(mod(seed-1,Nx*Ny)>=Ny);
    bw_x=sel_seed-Ny;
    ToUpdate=isnan(f_new(bw_x,1));
    new_seed=[new_seed; bw_x(ToUpdate)];
    f_new(bw_x(ToUpdate),:)=f_new(sel_seed(ToUpdate),:);

    sel_seed=seed(seed<=nel-Ny*Nx);
    fw_z=sel_seed+Ny*Nx;
    ToUpdate=isnan(f_new(fw_z,1));
    new_seed=[new_seed; fw_z(ToUpdate)];
    f_new(fw_z(ToUpdate),:)=f_new(sel_seed(ToUpdate),:);

    sel_seed=seed(seed>Ny*Nx);
    bw_z=sel_seed-Ny*Nx;
    ToUpdate=isnan(f_new(bw_z,1));
    new_seed=[new_seed; bw_z(ToUpdate)];
    f_new(bw_z(ToUpdate),:)=f_new(sel_seed(ToUpdate),:);

    seed=new_seed;

    it=it+1;

end

