function [ Fc] = projectF( F )
[V,D]=eig(F);

dd=diag(D);
[~,inds]=sort(abs(dd),'descend');
d=dd(inds(1:6));
pd=d(d>0);
nd=d(d<0);

VV=V(:,inds(1:6));
VP=VV(:,d>0);
VN=VV(:,d<0);
DDN=diag(nd);
DDP=diag(pd);


Fc=VP*DDP*VP'+VN*DDN*VN';
end

