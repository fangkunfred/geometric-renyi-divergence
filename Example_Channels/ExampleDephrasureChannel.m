% compute the Choi matrix of the dephrasure channel

function JN = ExampleDephrasureChannel(p,q)
I = eye(2);
JN = 0;
Z = [1 0; 0 -1];
for i = 1:2
    for j = 1:2
        rho = I(:,i)*I(:,j)';
        Noutput = blkdiag((1-q)*((1-p)*rho+p*Z*rho*Z), q*trace(rho));
        JN = JN + Tensor(rho,Noutput);
    end
end

end