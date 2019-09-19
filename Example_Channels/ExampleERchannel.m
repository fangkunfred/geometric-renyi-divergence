% Choi matrix of the quantum erasure channel with parameter p and dimension
% d
function J = ExampleERchannel(d,p)
 I1 = eye(d);
    I2 = eye(d+1);
    
    J1 = 0;
    for i = 1:d
        for j = 1:d
            J1 = J1 + Tensor(I1(:,i),I2(:,i))*Tensor(I1(:,j),I2(:,j))';
        end
    end
    J2 = 0;
    for i = 1:d
            J2 = J2 + Tensor(I1(:,i),I2(:,d+1))*Tensor(I1(:,i),I2(:,d+1))';
    end
    J = (1-p)*J1 + p*J2;
end