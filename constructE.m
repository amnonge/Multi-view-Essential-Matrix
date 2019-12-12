function [err,ourR,ourT,camList,newE] = constructE(Xs,Edges,Triplets)
Rs = {};
Ts = {};
for i = 1:length(Xs)
    Xs{i} = Xs{i}/norm(Xs{i});
    [R_output,T_output] = calculateTi(projectF( Xs{i} ) ,size(Xs{i},1)/3);
    Rs{i} = R_output;
    Ts{i} = T_output;
end
err = 0;
for i = 1:size(Edges,1)
    curEdge = Edges(i,:);
    [C,ia,ib] = intersect(Triplets(curEdge(1),1:3),Triplets(curEdge(2),1:3));
    [P,S,D] = findGlobalRotation(Rs{curEdge(1)}{ia(1)},Rs{curEdge(1)}{ia(2)},Rs{curEdge(2)}{ib(1)},Rs{curEdge(2)}{ib(2)},Ts{curEdge(1)}{ia(1)},Ts{curEdge(1)}{ia(2)},Ts{curEdge(2)}{ib(1)},Ts{curEdge(2)}{ib(2)});
    
    Rs{curEdge(2)}{1} = P*Rs{curEdge(2)}{1};
    Ts{curEdge(2)}{1} = S*P*Ts{curEdge(2)}{1}-D;
    
    Rs{curEdge(2)}{2} = P*Rs{curEdge(2)}{2};
    Ts{curEdge(2)}{2} = S*P*Ts{curEdge(2)}{2}-D;
    
    Rs{curEdge(2)}{3} = P*Rs{curEdge(2)}{3};
    Ts{curEdge(2)}{3} = S*P*Ts{curEdge(2)}{3}-D;
    
    err = err + norm(Ts{curEdge(1)}{ia(1)}-Ts{curEdge(2)}{ib(1)}) + norm(Ts{curEdge(1)}{ia(2)}-Ts{curEdge(2)}{ib(2)});
end

Us = [];
Vs = [];

ourR = {};
ourT = {};
camList = unique(Triplets);
for j = 1:length(camList) 
    i = camList(j);
    [is,js] = find(Triplets==i);
    Vs = [Vs;Rs{is(1)}{js(1)}'];
    Us = [Us;Rs{is(1)}{js(1)}'*getCrossM(Ts{is(1)}{js(1)})];
    
    ourR{end+1} = Rs{is(1)}{js(1)};
    ourT{end+1} = Ts{is(1)}{js(1)};
end
newE = Us*Vs'+Vs*Us';
end