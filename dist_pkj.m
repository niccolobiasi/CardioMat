function distM=dist_pkj(nodes,elem)
%distM=dist_pkj(nodes,elem) returns the distances between nodes in the
%purkinje tree defined by nodes and elements. DistM is a symmetric matrix
%NnodesxNnodes (Nnodes= number of nodes in the tree). Each element (i,j)
%represents the distance between node i and node j.
N_nodes=size(nodes,1);
distM=NaN(N_nodes);
for k=1:N_nodes
    curr_node=k;
    act_times=NaN(N_nodes,1);
    act_times(curr_node)=0;
    while nnz(isnan(act_times))
        [elem_in, loc]=ismember(elem,curr_node);
        loc=loc(elem_in);
        [seg,points]=ind2sub(size(elem_in),find(elem_in));
        uptCol=NaN(size(points));
        uptCol(points==1)=2;
        uptCol(points==2)=1;
        nodesToact=elem(sub2ind(size(elem),seg,uptCol));
        toRem=~isnan(act_times(nodesToact));
        nodesToact(toRem)=[];
        loc(toRem)=[];

        pathLen=vecnorm(nodes(nodesToact,:)-nodes(curr_node(loc),:),2,2);
        act_times(nodesToact)=act_times(curr_node(loc))+pathLen;

        curr_node=nodesToact;
    end


    distM(k,:)=act_times';
end

