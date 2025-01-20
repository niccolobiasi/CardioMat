function [aha, points]=ahaSegments(VoxelMat,apicobasal,rotational)
% aha=ahaSegments(VoxelMat,apicobasal,rotational) returns a
% matrix (aha) with same size of VoxelMat indicating AHA segments on the
% geometry described by VoxelMat. 
%[aha, points]=ahaSegments(VoxelMat,apicobasal,rotational) also returns the
%centroid of each segment.
% The function requires apicobasal and rotational coordinates.

aha=nan(nnz(VoxelMat),1);
aha(apicobasal==0)=17;
aha(apicobasal>0.66 & rotational>-2*pi/3 & rotational<=-pi/3)=1;
aha(apicobasal>0.66 & rotational>-pi/3 & rotational<=0)=2;
aha(apicobasal>0.66 & rotational>0 & rotational<=pi/3)=3;
aha(apicobasal>0.66 & rotational>pi/3 & rotational<=2*pi/3)=4;
aha(apicobasal>0.66 & rotational>2*pi/3 & rotational<=pi)=5;
aha(apicobasal>0.66 & rotational>=-pi & rotational<=-2*pi/3)=6;

aha(apicobasal>0.33 & apicobasal<=0.66  & rotational>-2*pi/3 & rotational<=-pi/3)=7;
aha(apicobasal>0.33 & apicobasal<=0.66 & rotational>-pi/3 & rotational<=0)=8;
aha(apicobasal>0.33 & apicobasal<=0.66 & rotational>0 & rotational<=pi/3)=9;
aha(apicobasal>0.33 & apicobasal<=0.66 & rotational>pi/3 & rotational<=2*pi/3)=10;
aha(apicobasal>0.33 & apicobasal<=0.66 & rotational>2*pi/3 & rotational<=pi)=11;
aha(apicobasal>0.33 & apicobasal<=0.66 & rotational>=-pi & rotational<=-2*pi/3)=12;

aha(apicobasal>0.0 & apicobasal<=0.33  & rotational>=-3*pi/4 & rotational<=-pi/4)=13;
aha(apicobasal>0.0 & apicobasal<=0.33  & rotational>-pi/4 & rotational<=pi/4)=14;
aha(apicobasal>0.0 & apicobasal<=0.33  & rotational>pi/4 & rotational<=3*pi/4)=15;
aha(apicobasal>0.0 & apicobasal<=0.33  & abs(rotational)>3*pi/4)=16;

if nargout>1
    points=nan(17,3);
    [ii, jj, kk]=ind2sub(size(VoxelMat),find(VoxelMat));

    apic_dist=(apicobasal-0.83).^2;

    [~,ind]=min(apic_dist+(rotational+pi/2).^2);
    points(1,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational+pi/6).^2);
    points(2,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational-pi/6).^2);
    points(3,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational-pi/2).^2);
    points(4,:)=[jj(ind) ii(ind) kk(ind)];
    
    [~,ind]=min(apic_dist+(rotational-5*pi/6).^2);
    points(5,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational+5*pi/6).^2);
    points(6,:)=[jj(ind) ii(ind) kk(ind)];

    apic_dist=(apicobasal-0.5).^2;

    [~,ind]=min(apic_dist+(rotational+pi/2).^2);
    points(7,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational+pi/6).^2);
    points(8,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational-pi/6).^2);
    points(9,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational-pi/2).^2);
    points(10,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational-5*pi/6).^2);
    points(11,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational+5*pi/6).^2);
    points(12,:)=[jj(ind) ii(ind) kk(ind)];

    apic_dist=(apicobasal-0.166).^2;

    [~,ind]=min(apic_dist+(rotational+pi/2).^2);
    points(13,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational).^2);
    points(14,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(rotational-pi/2).^2);
    points(15,:)=[jj(ind) ii(ind) kk(ind)];

    [~,ind]=min(apic_dist+(abs(rotational)-pi).^2);
    points(16,:)=[jj(ind) ii(ind) kk(ind)];

    
    points(17,:)=[mean(jj(aha==17)) mean(ii(aha==17)) mean(kk(aha==17))];
    
end
        