function [YS,Y,Xs1] = Optimization( FN,Cf,tts)
nn=size(FN,1);
n=size(Cf,1);

Ysize = max(max(Cf));

GAMMAS1=cell(n,1);
GAMMAS2=cell(n,1);
FNS=cell(n,1);
YS=cell(n,1);
for i=1:n
    GAMMAS1{i}=zeros(9,9);
    GAMMAS2{i}=zeros(9,9);
    FNS{i}=getTrippleInds(FN, Cf(i,:) );
    YS{i}=FNS{i};
end

[MapMatrix,sizeOfCell] = getMapMat(FNS,Cf,Ysize);

for i=1:100
    tempCell1 = cellfun(@substruction,YS,GAMMAS1,'UniformOutput', false);
    tempCell2 = cellfun(@substruction,YS,GAMMAS2,'UniformOutput', false);
    Xs1 = cellfun(@proxSVD,tempCell1,tempCell1,'UniformOutput', false);
    Xs2 = cellfun(@proxSVDNice,tempCell2,tts,'UniformOutput', false);
    Y=solve_Y4(GAMMAS1,Xs1,nn,Cf,FNS,sizeOfCell,MapMatrix,Xs2,GAMMAS2);
    
    for j=1:n
        YS{j}=genMatFromCell(Y,Cf(j,:),MapMatrix);
        GAMMAS1{j}=  (GAMMAS1{j}+Xs1{j}-YS{j});
        GAMMAS2{j}=  (GAMMAS2{j}+Xs2{j}-YS{j});
    end
end
end





function Y=solve_Y4(GAMMAS1,Xs1,nn,Cf,FNS,sizeOfCell,MapMatrix,Xs2,GAMMAS2)
Y=cell(1,sizeOfCell);
W=cell(1,sizeOfCell);
param1 = 0.1;
param2 = 100;
for i=1:length(GAMMAS1)
    Y = insertData(Y,Cf(i,:),(Xs1{i}+GAMMAS1{i})/param1+(Xs2{i}+GAMMAS2{i})/param2+FNS{i},MapMatrix);
    W = insertDataW(W,Cf(i,:),(1+1/param1 + 1/param2),MapMatrix);
end
Y = cellfun(@division,Y,W,'UniformOutput',false);
Y = cellfun(@SVDEsmall,Y,'UniformOutput',false);

end







function F=proxSVD( Y1, Y2 )

res=10*Y1/11+Y2/11;
[W,D] = eigs(0.5*(res+res'));
[~,c1] = sort(diag(D),'descend');
X = W(:,c1(1:3));
[~,c] = sort(diag(D),'ascend');
Y = W(:,c(1:3));
D(c(1:3),c(1:3)) = (D(c(1:3),c(1:3)) -D(c1(1:3),c1(1:3)))*0.5;
V = sqrt(0.5)*( X +  Y);
U = sqrt(0.5)*( X -  Y);
F= V*(U*(-D(c(1:3),c(1:3))))'+(U*(-D(c(1:3),c(1:3))))*V';

end



function [ma,tr]=getTrippleInds(BigMatrix, tripleInds )
tr=[tripleInds(1)*3-2:tripleInds(1)*3 tripleInds(2)*3-2:tripleInds(2)*3 tripleInds(3)*3-2:tripleInds(3)*3];
ma=BigMatrix(tr,tr);
end

function tr=getTrippleInds2(tripleInds )
tr=[tripleInds(1)*3-2:tripleInds(1)*3 tripleInds(2)*3-2:tripleInds(2)*3 tripleInds(3)*3-2:tripleInds(3)*3];
end

function Esmall = SVDEsmall(Esmall)

if ~isempty(Esmall)
    [u,d,v]=svd(Esmall);
    ee=(d(2,2)+d(1,1))/2;
    d(1,1)=ee;
    d(2,2)=ee;
    Esmall=u*d*v';
end

end


function tripMat = genMatFromCell(BigCell,trip,MapMatrix)

tripMat = zeros(9);
tripMat(1:3,4:6) = BigCell{MapMatrix(trip(1),trip(2))};
tripMat(1:3,7:9) = BigCell{MapMatrix(trip(1),trip(3))};
tripMat(4:6,7:9) = BigCell{MapMatrix(trip(2),trip(3))};

tripMat = tripMat+tripMat';
end


function Y = insertData(Y,trip,data,MapMatrix)


if ~isempty(Y{MapMatrix(trip(1),trip(2))})
    Y{MapMatrix(trip(1),trip(2))} = Y{MapMatrix(trip(1),trip(2))}+data(1:3,4:6);
else
    Y{MapMatrix(trip(1),trip(2))}=data(1:3,4:6);
end

if ~isempty(Y{MapMatrix(trip(1),trip(3))})
    Y{MapMatrix(trip(1),trip(3))} = Y{MapMatrix(trip(1),trip(3))}+data(1:3,7:9);
else
    Y{MapMatrix(trip(1),trip(3))} = data(1:3,7:9);
end


if ~isempty(Y{MapMatrix(trip(2),trip(3))})
    Y{MapMatrix(trip(2),trip(3))} = Y{MapMatrix(trip(2),trip(3))}+data(4:6,7:9);
else
    Y{MapMatrix(trip(2),trip(3))} = data(4:6,7:9);
end

end

function Y = insertDataW(Y,trip,data,MapMatrix)

if ~isempty(Y{MapMatrix(trip(1),trip(2))})
    Y{MapMatrix(trip(1),trip(2))} = Y{MapMatrix(trip(1),trip(2))}+data;
else
    Y{MapMatrix(trip(1),trip(2))}=data;
end

if ~isempty(Y{MapMatrix(trip(1),trip(3))})
    Y{MapMatrix(trip(1),trip(3))} = Y{MapMatrix(trip(1),trip(3))}+data;
else
    Y{MapMatrix(trip(1),trip(3))} = data;
end


if ~isempty(Y{MapMatrix(trip(2),trip(3))})
    Y{MapMatrix(trip(2),trip(3))} = Y{MapMatrix(trip(2),trip(3))}+data;
else
    Y{MapMatrix(trip(2),trip(3))} = data;
end

end


function res = division(X,Y)
if ~isempty(Y)
    res = X./Y;
else
    res = X;
end
end

function res = substruction(X,Y)
res = X-Y;
end

function [MapMatrix,counter] = getMapMat(FNs,trips,Ysize)
tempBigMat = zeros(Ysize);
for i=1:length(FNs)
    tr = getTrippleInds2(trips(i,:) );
    tempBigMat(tr,tr)=FNs{i};
end

tempBigMat = tempBigMat~=0;
tempBigMat = tempBigMat(1:3:end,1:3:end);
tempBigMat = tril(tempBigMat);

counter = 1;
mapBigMat = zeros(size(tempBigMat));
for i = 1:size(tempBigMat,1)*size(tempBigMat,2)
    if tempBigMat(i)
        mapBigMat(i) = counter;
        counter = counter+1;
    end
end
MapMatrix = mapBigMat+mapBigMat';



end

function F=proxSVDNice( Y1,ttsVal )

res = Y1;

if ttsVal<0
    [u,d,v]=svd(res);
    temp=diag(d);
    temp(7:end)=0;
    tempp1=0.5*(temp(1)+temp(2));
    
    tempp2=0.5*(temp(3)+temp(4));
    
    tempp3=0.5*(temp(5)+temp(6));
    temp(1)=tempp1;
    temp(2)=tempp1;
    temp(3)=tempp2;
    temp(4)=tempp2;
    temp(5)=tempp3;
    temp(6)=tempp3;
    
    d=diag(temp);
    F = u*d*v';
    
else
    
    for kkk = 1:1
        res = eigenF(0.5*(res+res'));
        
    end
    
    F = res;
    
end

end


function F = eigenF(res)

[W,D] = eigs(res);
[~,c1] = sort(diag(D),'descend');
X = W(:,c1(1:3));
[~,c] = sort(diag(D),'ascend');
Y = W(:,c(1:3));
minVal = -inf;
signs = [-1,1];


for i1 = signs
    for i2 = signs
        for i3 = signs
            tempX = X*diag([i1,i2,i3]);
            tempY = Y ;
            tempV = 0.5*(tempX + tempY);
            V1 = tempV(1:3,1:3)'*tempV(1:3,1:3);
            V2 = tempV(4:6,1:3)'*tempV(4:6,1:3);
            V3 = tempV(7:9,1:3)'*tempV(7:9,1:3);
            if (norm(diag(V1))/norm(V1,'fro')+norm(diag(V2))/norm(V2,'fro')+norm(diag(V3))/norm(V3,'fro'))>minVal         %abs(V12)<treshold & abs(V13)<treshold & abs(V23)<treshold
                minVal =(norm(diag(V1))/norm(V1,'fro')+norm(diag(V2))/norm(V2,'fro')+norm(diag(V3))/norm(V3,'fro'));
                V = sqrt(0.5)*(tempX + tempY);
                U = sqrt(0.5)*(tempX - tempY);
            end
        end
    end
end

Vn=V;

[uu,dd,vv]=svd(V(1:3,1:3));
Vn(1:3,1:3)=uu*diag(repmat(mean(diag(dd)),3,1))*vv';
[uu,dd,vv]=svd(V(4:6,1:3));
Vn(4:6,1:3)=uu*diag(repmat(mean(diag(dd)),3,1))*vv';
[uu,dd,vv]=svd(V(7:9,1:3));
Vn(7:9,1:3)=uu*diag(repmat(mean(diag(dd)),3,1))*vv';
X = (Vn+U)*sqrt(0.5);
Y =  (Vn-U)*sqrt(0.5);
F = X*D(c1(1:3),c1(1:3))*X'+Y*D(c(1:3),c(1:3))*Y';




if min(norm(F/norm(F)-res/norm(res)),norm(F/norm(F)+res/norm(res)))> 0.5%0.005
    
    [u,d,v]=svd(res);
    temp=diag(d);
    temp(7:end)=0;
    tempp1=0.5*(temp(1)+temp(2));
    
    tempp2=0.5*(temp(3)+temp(4));
    
    tempp3=0.5*(temp(5)+temp(6));
    temp(1)=tempp1;
    temp(2)=tempp1;
    temp(3)=tempp2;
    temp(4)=tempp2;
    temp(5)=tempp3;
    temp(6)=tempp3;
    
    d=diag(temp);
    F = u*d*v';
    
end
end

