function [ firstGroup,tempGraph] = reduceTriplets(c,score,G)
names = cell(1,size(c,1));
for i = 1:size(c,1)
    names{i} = num2str(i);
end
G.Nodes.Name = names';
[~,idds] = sort(score);

AA=adjacency(G);
[is,js]=find(AA);
keep=(js>is);
EEE=[is,js];
EEE=EEE(keep,:);

EE=EEE';
edges=EE(:);
cc=c';
[allCams,~,tr]=unique(cc(:));
numTriplets=size(c,1);
numcam=length(allCams);
mapTriangleIndices=0:numTriplets-1;
mapTriangleIndices2=0:numTriplets-1;
mask=ones(1,numTriplets);
mask2=ones(1,numTriplets);
mask3=ones(1,numTriplets);
numedges=size(EE,2);

out2=formexgraph6(int32(idds),int32(tr),int32(allCams),int32(mapTriangleIndices),int32(edges),numTriplets,numcam,numedges,logical(mask),int32(mapTriangleIndices2),logical(mask2),logical(mask3));
firstGroup = out2>0;
tempGraph = subgraph(G,find(firstGroup));
end

