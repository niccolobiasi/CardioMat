function phi= computePhase(VoxelMat,dx,dt,xi,useGPU)
% phi= computePhase(VoxelMat,dx,dt,xi,useGPU) returns the phase field phi computed for the
% voxelized geometry VoxelMat. dx, dt, and xi are mandatory parameters
% used for the computation of the phase field (see Biasi et al, Plos one
% 2023 for details). useGPU specifies if the function should run on GPU
% (true) or on CPU (false).
if xi==0
    phi=double(VoxelMat);
    return
end

dx2=dx^2;
T=20;
N=T/dt+1;

[FV,extInd]=computeSurface(VoxelMat,dx);

%diffuse state + thresholding
phi=imgaussfilt3(double(VoxelMat),8*xi/dx);
mask=phi>(mean(phi,'all')*0.5);
phi=zeros(size(VoxelMat));
phi(VoxelMat & mask)=1;
%laplace matrix specific
S=laplaceMatrix(mask);

xi2=xi^2;
if useGPU
    dt=gpuArray(dt);
    dx2=gpuArray(dx2);
    xi2=gpuArray(xi2);
    phi_int=gpuArray(phi(mask));
    S=gpuArray(S);

    for i=1:N

                zz=2*phi_int-1;
                phi_int=phi_int+dt*(xi2/dx2*(S*phi_int)-(zz.^3-zz));
    end
    phi(mask) = gather(phi_int);
else
    for i=1:N

        zz=2*phi-1;
        phi=phi+dt*(xi2/dx2*(S*phi)-(zz.^3-zz));

    end
end


phi=reshape(phi,size(VoxelMat,1),size(VoxelMat,2),size(VoxelMat,3));
end
