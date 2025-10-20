function VoxelMat=cleanGeometry(VoxelMat)
% VoxelMat=cleanGeometry(VoxelMat) remove unreachable voxels from the
% input geometry VoxelMat. Only the largest domain is mantained.
eik1=-1*ones(numel(VoxelMat),1);
k=1;
domain=zeros(size(VoxelMat));
bs=0;
while nnz(eik1(VoxelMat)<0)     
    seed=find(eik1==-1 & VoxelMat(:));
    if length(seed)<bs
        break;
    end
    seed=seed(randi(length(seed)));
    eik1=solveEikonal(VoxelMat,seed);
    sel=eik1>=0;
    Ns=nnz(sel);
    if Ns>bs
        bk=k;
        bs=Ns;
    end
    domain(sel)=k;
    k=k+1;
end
VoxelMat(domain~=bk)=0;

