function phi= computePhase(VoxelMat,dx,dt,xi,useGPU)
% computePhase returns the phase field computed for the
% voxelized geometry VoxelMat. dx, dt, and xi are mandatory parameters
% used for the computation of the phase field (see Biasi et al, Plos one
% 2023 for details). useGPU specifies if the function should run on GPU
% (true) or on CPU (false).

phi=zeros(size(VoxelMat));
phi(VoxelMat)=1;
phi=phi(:);
dx2=dx^2;
T=20;
N=T/dt+1;
S=laplaceMatrix(size(VoxelMat));

xi2=xi^2;
if useGPU
    dt=gpuArray(dt);
    dx2=gpuArray(dx2);
    xi2=gpuArray(xi2);
    phi=gpuArray(phi);
        S=gpuArray(S);

    for i=1:N

                zz=2*phi-1;
                phi=phi+dt*(xi2/dx2*(S*phi)-(zz.^3-zz));

    end
    phi = gather(phi);
else
    for i=1:N

        zz=2*phi-1;
        phi=phi+dt*(xi2/dx2*(S*phi)-(zz.^3-zz));

    end
end


phi=reshape(phi,size(VoxelMat,1),size(VoxelMat,2),size(VoxelMat,3));
end
