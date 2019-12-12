function [RF1,RF2,T1,T2,RD1,RD2,c_opt2]= getBestError(groundRs,Rs,groundTs,ts,indices)

[Rmat,Tmat] = convertListtoMat(Rs,ts);
[RmatGT,TmatGT] = convertListtoMat(groundRs,groundTs);


[SO3Mats_corr2, MSE_rots2, R_global] = GlobalSOdCorrectLeft(Rmat, RmatGT);
[t_fit2,t_opt2,c_opt2,NRMSE_LUD2,~] = SimpleTransScaleRemove(R_global*Tmat, TmatGT,'L1');


%% Translation Error
T1=mean(sqrt(sum((t_fit2 - TmatGT).^2)));
T2=median(sqrt(sum((t_fit2 - TmatGT).^2)));


%% Rotation Error in degrees
% Delete those cameras with f < 10^-5 for comparison

[E_res,e,normE]=CompareRotations(SO3Mats_corr2 ,RmatGT);
RD1=E_res(1);
RD2=E_res(2);
RF1=normE(1);
RF2=normE(2);

end

function [Rmat,Tmat]=convertListtoMat(RList,TList)

Rmat = zeros(3,3,length(RList));
Tmat = zeros(3,length(TList));

for i = 1:length(RList)
    Rmat(:,:,i) = RList{i};
    Tmat(:,i) = TList{i};
end
    

end
