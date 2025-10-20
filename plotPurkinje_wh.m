function plotPurkinje_wh(ventr,atr,transmural,phi,pkn_nodes,pkn_elem,pkn_times)

if length(transmural)==numel(ventr)
    transmural =transmural(ventr);
end

if length(phi)==numel(atr)
    phi =phi(atr);
end

VoxelMat=atr | ventr;
phi_wh=nan(numel(VoxelMat),1);
phi_wh(atr(:))=phi;
phi_wh(ventr(:))=transmural;

plotPurkinje(VoxelMat,phi_wh,pkn_nodes,pkn_elem,pkn_times)