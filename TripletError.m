function [ error ] = TripletError(curRR )
R12=curRR(1:3,4:6);
R13=curRR(1:3,7:9);
R23=curRR(4:6,7:9);
error = norm(R12*R23*(R13')-eye(3))+1;
end

