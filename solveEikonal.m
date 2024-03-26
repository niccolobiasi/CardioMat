function eik=solveEikonal(VoxelMat, seed)
% solveEikonal solves the Eikonal equation by simulating wavefront
% propagation in the domain. 
% eik=solveEikonal(VoxelMat, seed) outputs the steps (in terms of number of
% voxels) needed to achieve the point inside the voxelized geometry
% starting from the seed. If there are unreacheable voxels the function
% notifies it to the user
%Note: the function assumes VoxelMat is equal to 0 in the computational
%border
%
eik=-1*ones(numel(VoxelMat),1);
eik(seed)=0;
nin=nnz(VoxelMat);
[Ny,Nx,~]=size(VoxelMat);

while length(find(eik>-1))<nin

    if isempty(seed)
        break;
    end
    new_seed=[];



    fw_y=seed+1;
    ToUpdate=VoxelMat(fw_y) & eik(fw_y)==-1;
    new_seed=[new_seed; fw_y(ToUpdate)];
    eik(fw_y(ToUpdate))=eik(seed(ToUpdate))+1;

    bw_y=seed-1;
    ToUpdate=VoxelMat(bw_y) & eik(bw_y)==-1;
    new_seed=[new_seed; bw_y(ToUpdate)];
    eik(bw_y(ToUpdate))=eik(seed(ToUpdate))+1;

    fw_x=seed+Ny;
    ToUpdate=VoxelMat(fw_x) & eik(fw_x)==-1;
    new_seed=[new_seed; fw_x(ToUpdate)];
    eik(fw_x(ToUpdate))=eik(seed(ToUpdate))+1;

    bw_x=seed-Ny;
    ToUpdate=VoxelMat(bw_x) & eik(bw_x)==-1;
    new_seed=[new_seed; bw_x(ToUpdate)];
    eik(bw_x(ToUpdate))=eik(seed(ToUpdate))+1;

    fw_z=seed+Ny*Nx;
    ToUpdate=VoxelMat(fw_z) & eik(fw_z)==-1;
    new_seed=[new_seed; fw_z(ToUpdate)];
    eik(fw_z(ToUpdate))=eik(seed(ToUpdate))+1;

    bw_z=seed-Ny*Nx;
    ToUpdate=VoxelMat(bw_z) & eik(bw_z)==-1;
    new_seed=[new_seed; bw_z(ToUpdate)];
    eik(bw_z(ToUpdate))=eik(seed(ToUpdate))+1;


    seed=new_seed;
end