function [pkn_nodes,pkn_elem,pkn_times]=insertAVN(pkn_nodes,pkn_elem,atr,pkn_times,res,cond_vel)
% [pkn_nodes,pkn_elem,pkn_times]=insertAVN(pkn_nodes,pkn_elem,atr,pkn_times)
% inserts a user selected AVN point at the beginning of the Purkinje tree.
% The AVN point is selected from the atrial surface.
% The His starting point in the original network is substituted with the
% mean point across the AVN node, the left ventricle network starting point
%, and the right ventricle network starting point.
% The function assumes the order of the elements in pkn_elem is the same as
% produced by createPurkinje_biv.m.
% Optional inputs are the resolution (res) and the conduction velocity
% (cond_vel, needed fr adjusting activation time only).

if nargin<6
    cond_vel=340;
end

if nargin<5
    res=0.25;
elseif isempty(res)
    res=0.25;
end

[FV_atr,extInd_atr]=computeSurface(atr,res);
[~, avn_ind]=selectPoint(atr,FV_atr,extInd_atr,res);
[avn_y, avn_x, avn_z]=ind2sub(size(atr),avn_ind);
pkn_elem=[1 2; pkn_elem+1];
pkn_nodes=[avn_x avn_y avn_z; pkn_nodes];
pkn_nodes(2,:)=(pkn_nodes(1,:)+pkn_nodes(3,:)+pkn_nodes(pkn_elem(3,2),:))/3;

if nargout>2
    if nargin<4 
        error('Please provide pkn_times input');
    elseif isempty(pkn_times)
        error('Please provide pkn_times input');
    end
     dist1=norm(pkn_nodes(2,:)-pkn_nodes(1,:));
     distL=norm(pkn_nodes(2,:)-pkn_nodes(3,:));
     distR=norm(pkn_nodes(2,:)-pkn_nodes(pkn_elem(3,2),:));
     pkn_times(2:(pkn_elem(3,2)-2))=pkn_times(2:(pkn_elem(3,2)-2))-pkn_times(2)+distL*res/cond_vel*100;
     pkn_times((pkn_elem(3,2)-1):end)=pkn_times((pkn_elem(3,2)-1):end)-pkn_times(pkn_elem(3,2)-1)+distR*res/cond_vel*100;
     pkn_times=[0; pkn_times+dist1*res/cond_vel*100];
end

end