function [tripList,Vlist,conCompList,ttsList] = buildTriplets(pointMatchesInliers,FN,Hmat,TijMat,ParamOpt)

numcam=size(pointMatchesInliers,2);



inliersNum= pointMatchesInliers;
G=getViewingGraph(inliersNum,  numcam,FN);
tripletViewingGraph=getTripletViewingGraph(G);
[c]=extractTripletsFromViewingGraph(graph(tripletViewingGraph));
[transCons,tts,Inliersall,tripletsErrors]=characteriseTriplets(TijMat,FN,c,Hmat,inliersNum);
tts=1./tts;


tripletsErrors = 1./tripletsErrors;

extra = 0;

colTresh = 1/ParamOpt.col;
colTresh2 = 1/(ParamOpt.rot+extra);
colTresh3 = ParamOpt.tran;

maskVec = tts>colTresh & tripletsErrors(:,1)>colTresh2 & transCons<colTresh3; % & tripletsErrors(:,2)>colTresh3; %& tripletsErrors(:,1)>colTresh2;%& tripletsErrors(:,2)>colTresh3;
c=c(maskVec ,:);
tripletsErrors=tripletsErrors(maskVec ,:);
tts=tts(maskVec,:);


Gt=getTripletGraph(c);



conCompArray = conncomp(Gt);
tts2=  tripletsErrors(:,1);

unList = unique(conCompArray);
argmax = 0;
indmax = 0;
for i = 1:length(unList)
    temp = length(unique(c(conCompArray==unList(i))));
    if temp>argmax
        argmax = temp;
        indmax = unList(i);
    end
end



tempidx = find(indmax==conCompArray);
tempGraph = subgraph(Gt,tempidx);
[ firstGroup,tempGraph] = reduceTriplets(c(tempidx,:),tts2(tempidx,:),tempGraph);
tempc=c(tempidx,:);
tempc = tempc(firstGroup,:);
temptts = tts(tempidx);
temptts = temptts(firstGroup);

tempV = bfsearch(tempGraph,1,{'edgetonew'});
ttsList = temptts;
conCompList = tempGraph;
Vlist = tempV;
tripList = tempc;


end




function G=getViewingGraph(inliersNum,  numcam,FN)
curF = abs(FN)>0;
curF = curF(1:3:end,1:3:end);
inliersNum = inliersNum+inliersNum';
adjec =  inliersNum>15 & curF;
adjec = adjec./inliersNum;
adjec(isnan(adjec))=0;
G=graph(adjec);
end

function [transCons,tts,Inliersall,tripletsErrors]=characteriseTriplets(TijMat,FN,c,Hmat,inliersNum)
transCons=zeros(size(c,1),1);
tts=zeros(size(c,1),1);
Inliersall=zeros(size(c,1),3);
tripletsErrors=zeros(size(c,1),2);
for i=1:size(c,1)
    try
        tts(i)= max([abs(TijMat{c(i,1),c(i,2)}'* TijMat{c(i,1),c(i,3)}/(norm(TijMat{c(i,1),c(i,2)})*norm(TijMat{c(i,1),c(i,3)}))),abs(TijMat{c(i,2),c(i,1)}'* TijMat{c(i,2),c(i,3)}/(norm(TijMat{c(i,2),c(i,1)})*norm(TijMat{c(i,2),c(i,3)}))),abs(TijMat{c(i,3),c(i,1)}'* TijMat{c(i,3),c(i,2)}/(norm(TijMat{c(i,3),c(i,1)})*norm(TijMat{c(i,3),c(i,2)})))]);
        transCons(i) = abs(acos(TijMat{c(i,1),c(i,2)}'* TijMat{c(i,1),c(i,3)}/(norm(TijMat{c(i,1),c(i,2)})*norm(TijMat{c(i,1),c(i,3)})))+acos(TijMat{c(i,2),c(i,1)}'* TijMat{c(i,2),c(i,3)}/(norm(TijMat{c(i,2),c(i,1)})*norm(TijMat{c(i,2),c(i,3)})))+acos(TijMat{c(i,3),c(i,1)}'* TijMat{c(i,3),c(i,2)}/(norm(TijMat{c(i,3),c(i,1)})*norm(TijMat{c(i,3),c(i,2)})))-pi);
    catch
        tts(i) = 1;
        transCons(i) = 1;
    end
    inds=[3*c(i,1)-2:3*c(i,1) 3*c(i,2)-2:3*c(i,2) 3*c(i,3)-2:3*c(i,3) ];
    curFF=FN(inds,inds);
    curRR = Hmat(inds,inds);
    [ curFF ] = normalizeForbineusNorm( curFF );
    
    [err]=TripletError(curRR );
    tripletsErrors(i,1)=err;

    Inliersall(i,1)=inliersNum(c(i,1),c(i,2));
    Inliersall(i,2)=inliersNum(c(i,1),c(i,3));
    Inliersall(i,3)=inliersNum(c(i,2),c(i,3));
end
end

function tripletGraph=getTripletViewingGraph(G)
T=minspantree(G);
tripletGraph=adjacency(T);
graphsTrips=cell(4,1);
graphsTrips{1}=graph(tripletGraph);
try
    for i=1:25
        G = rmedge(G,T.Edges.EndNodes(:,1),T.Edges.EndNodes(:,2));
        T=minspantree(G,'root',G.Edges.EndNodes(1));
        tripletGraph=tripletGraph+adjacency(T);
        graphsTrips{i+1}=graph(tripletGraph);
    end
catch
    'error occured in tree'
    i
end
end

function Gt=getTripletGraph(c)
adjec=buildTripletGraph(int32(c'));
Gt=graph(adjec);
end