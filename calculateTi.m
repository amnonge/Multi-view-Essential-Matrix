function [R_output,T_output] = calculateTi(E_est,N)
[W,D] = eigs(E_est);

[~,c] = sort(diag(D),'descend');
X = W(:,c(1:3));
[~,c] = sort(diag(D),'ascend');
Y = W(:,c(1:3));

minVal = 1;
signs = [-1,1];
for i1 = signs
    for i2 = signs
        for i3 = signs
            for i4 = signs
                for i5 = signs
                    for i6 = signs
                        tempX = X*diag([i1,i2,i3]);
                        tempY = Y*diag([i4,i5,i6]);
                        tempV = 0.5*(tempX + tempY);
                        
                        V12 = abs(dot(tempV(1:3,1),tempV(1:3,2)))+abs(dot(tempV(4:6,1),tempV(4:6,2)))+abs(dot(tempV(7:9,1),tempV(7:9,2)));
                        V13 = abs(dot(tempV(1:3,1),tempV(1:3,3)))+abs(dot(tempV(4:6,1),tempV(4:6,3)))+abs(dot(tempV(7:9,1),tempV(7:9,3)));
                        V23 = abs(dot(tempV(1:3,3),tempV(1:3,2)))+abs(dot(tempV(4:6,3),tempV(4:6,2)))+abs(dot(tempV(7:9,3),tempV(7:9,2)));
                        if max(abs([V12,V13,V23]))<minVal         %abs(V12)<treshold & abs(V13)<treshold & abs(V23)<treshold
                            minVal = max(abs([V12,V13,V23]));
                            V = 0.5*(tempX + tempY);
                            U = 0.5*(tempX - tempY);
                        end
                    end
                end
            end
        end
    end
end
alphavec = [];
for i = 1:N
    alphavec = [alphavec;nthroot(real(det(V(3*i-2:3*i,1:3))),3)];
end
ScalesMat = zeros(size(E_est));
for i = 1:length(alphavec)
    ScalesMat(3*i-2:3*i,3*i-2:3*i) = eye(3)*1/abs(alphavec(i));
end

E_est = ScalesMat*E_est*ScalesMat';

[u,s,v] = svd(E_est);

for i = 1:2:6
    entry = (s(i,i)+s(i+1,i+1))/2;
    s(i,i) = entry;
    s(i+1,i+1) = entry;
end

E_est = u*s*v';

[W,D] = eigs(E_est);

[~,c] = sort(diag(D),'descend');
X = W(:,c(1:3));
[~,c] = sort(diag(D),'ascend');
Y = W(:,c(1:3));

minVal = 1;
signs = [-1,1];
for i1 = signs
    for i2 = signs
        for i3 = signs
            for i4 = signs
                for i5 = signs
                    for i6 = signs
                        tempX = X*diag([i1,i2,i3]);
                        tempY = Y*diag([i4,i5,i6]);
                        tempV = 0.5*(tempX + tempY);
                        V12 = abs(dot(tempV(1:3,1),tempV(1:3,2)))+abs(dot(tempV(4:6,1),tempV(4:6,2)))+abs(dot(tempV(7:9,1),tempV(7:9,2)));
                        V13 = abs(dot(tempV(1:3,1),tempV(1:3,3)))+abs(dot(tempV(4:6,1),tempV(4:6,3)))+abs(dot(tempV(7:9,1),tempV(7:9,3)));
                        V23 = abs(dot(tempV(1:3,3),tempV(1:3,2)))+abs(dot(tempV(4:6,3),tempV(4:6,2)))+abs(dot(tempV(7:9,3),tempV(7:9,2)));
                        if max(abs([V12,V13,V23]))<minVal
                            minVal = max(abs([V12,V13,V23]));
                            V = (tempX + tempY);
                            U = 0.5*(tempX - tempY);
                        end
                    end
                end
            end
        end
    end
end

for i = 1:N
    
    [u,s,v] = svd(V(3*i-2:3*i,1:3));
end
R_output = cell(N,1);
T_output = cell(N,1);
for i = 1:N
    U(3*i-2:3*i,1:3) = inv(V(3*i-2:3*i,1:3))*U(3*i-2:3*i,1:3)*abs(D(c(1:3),c(1:3)));
    U(3*i-2:3*i,1:3) = 0.5*(U(3*i-2:3*i,1:3)-U(3*i-2:3*i,1:3)');
    tempT_i = U(3*i-2:3*i,1:3);
    [u,s,v] = svd(V(3*i-2:3*i,1:3));
    signOfV = sign(det(u*v'));
    V(3*i-2:3*i,1:3) = signOfV*u*v';
    signOfV = sign(det(V(3*i-2:3*i,1:3)));
    R_output{i} = signOfV*V(3*i-2:3*i,1:3)';
    T_output{i} = [tempT_i(3,2), tempT_i(1,3), tempT_i(2,1)]';
end


end

