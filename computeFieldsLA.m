function fields_LA=computeFieldsLA(LA,tag,res)
% fields_LA=computeFieldsLA(LA,tag,res) returns a cell array of 4 elements containing the
% distance fields needed for the assignement in the fiber orientation in
% the left atrium (voxelized geometry LA). 
% For further information see the related paper by Piersanti et al. (2021)
% The function requires the surface tags tag as produced by tagAtria
% function. The resolution is optional (default: 0.25 mm). Manual selection
% of the left atrial appendage is required.

if nargin<3 
    res=0.25;
end

[FV_l,extInd_l]=computeSurface(LA,res);
ind_surf_l=unique(extInd_l);
isSurf_l=false(size(LA));
isSurf_l(ind_surf_l)=1;
eik=solveEikonal(LA,find((tag==9 | (tag==0 & isSurf_l(:))) & LA(:)));
eik1=solveEikonal(LA,find(tag==8));
phi_l=eik1./(eik+eik1);
eik_lpv=solveEikonal(LA,find(tag==4));
eik_rpv=solveEikonal(LA,find(tag==5));
disp('Select left atrial appendage...')
la_app=selectSurface(LA,FV_l,extInd_l,res,1, true);
disp('Done!')
eik_lapp=solveEikonal(LA,find(la_app));
eik_mv=solveEikonal(LA,find(tag==6));

eik_lpv_rpv_mv=solveEikonal(LA,find(tag==4 | tag==5 | tag==6));
eik_lpv_rpv_lapp=solveEikonal(LA,find(tag==4 | tag==5 | la_app(:)));
eik_lpv_mv_lapp=solveEikonal(LA,find(tag==4 | tag==6 | la_app(:)));
eik_rpv_mv_lapp=solveEikonal(LA,find(tag==5 | tag==6 | la_app(:)));

psi_ab_la=(3*(eik_lapp./(eik_lapp+eik_lpv_rpv_mv))+...
    1*(eik_lpv./(eik_lpv+eik_rpv_mv_lapp))-1*eik_mv./(eik_mv+eik_lpv_rpv_lapp)-3*eik_rpv./(eik_rpv+eik_lpv_mv_lapp)+1)/2;

psi_v_la=eik_lpv./(eik_lpv+eik_rpv);

psi_r_la=eik_lpv_rpv_lapp./(eik_mv+eik_lpv_rpv_lapp);
fields_LA={phi_l,psi_ab_la,psi_v_la, psi_r_la};

end