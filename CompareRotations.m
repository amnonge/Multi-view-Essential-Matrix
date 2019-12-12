function [E,e,normE,norme]=CompareRotations(R1,R2)

e=zeros(size(R1,3),1); norme=zeros(size(R1,3),1);
for i=1:size(R1,3)
    if(sum(sum((R1(:,:,i)==0)))==9||any(any(isnan(R1(:,:,i))))||...
        sum(sum((R2(:,:,i)==0)))==9||any(any(isnan(R2(:,:,i)))));
        e(i,1)=NaN;
        norme(i,1)=NaN;
    else
        e(i,1)=acos(max(min((R1(1,:,i)*R2(1,:,i)'+R1(2,:,i)*R2(2,:,i)'+R1(3,:,i)*R2(3,:,i)'-1)/2,1),-1));
        norme(i,1)=norm(R1(:,:,i)-R2(:,:,i),'fro');
    end
end
e=e*180/pi;
i=~isnan(e);
E=[mean(e(i))  median(e(i)) sqrt(e(i)'*e(i)/sum(i))];
normE=[nanmean(norme) nanmedian(norme)];

end
