function [FN,M,R_gt,T_gt,K_gt] = getSubGraph(M,FN,R_gt,T_gt,K_gt,indices)

R_gt = R_gt(indices);
T_gt = T_gt(indices);
K_gt = K_gt(indices);

indicesx=indices*2-1;
indicesy=indices*2;
Mc=zeros(length(indices)*2,size(M,2));
Mc(1:2:end,:)=M(indicesx,:);
Mc(2:2:end,:)=M(indicesy,:);
M=Mc;

M=M(:,(sum(abs(M)>10^-5,1)>=4));
FNc=zeros(3*length(indices),3*length(indices));
for i=1:length(indices)-1
    for j=i+1:length(indices)
        FNc(3*i-2:3*i,3*j-2:3*j)=FN(3*indices(i)-2:3*indices(i),3*indices(j)-2:3*indices(j));
        FNc(3*j-2:3*j,3*i-2:3*i)=FNc(3*i-2:3*i,3*j-2:3*j)';
    end
end
FN=FNc;


end