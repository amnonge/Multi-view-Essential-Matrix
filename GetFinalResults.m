function [RF1,RF2,T1,T2,RD1,RD2,N]=GetFinalResults(dataset)
load(sprintf('Final_names/%s.mat',dataset));
load(sprintf('results/%s.mat',dataset));

A=cellfun(@(x) x(8:end),namesList,'UniformOutput',false);

[~,leftA,leftB]=intersect(A,files);
ourR=ourR(leftA);
ourT=ourT(leftA);

R_gt=R_gt(leftA);
T_gt=T_gt(leftA);


 [RF1,RF2,T1,T2,RD1,RD2,scale]= getBestError(R_gt,ourR,T_gt,ourT ,1:length(R_gt));

N=length(R_gt);


