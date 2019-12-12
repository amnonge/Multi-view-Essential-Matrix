function [P,S,D] = findGlobalRotation(Rs11,Rs12,Rs21,Rs22,Ts11,Ts12,Ts21,Ts22)
P = Rs11*Rs21';
S=(Ts11-Ts12)./(P*(Ts21-Ts22));
S=S(3);
Ts21s=S*(P*Ts21);
D=Ts21s-Ts11;
end