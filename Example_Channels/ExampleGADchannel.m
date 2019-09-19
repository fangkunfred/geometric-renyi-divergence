% Choi matrix of the generalized amplitude damping channel

function JN = ExampleGADchannel(r,N)
    A{1} = sqrt(1-N)*[1 0; 0 sqrt(1-r)];
    A{2} = sqrt(r*(1-N))*[0 1; 0 0];
    A{3} = sqrt(N)*[sqrt(1-r) 0; 0 1];
    A{4} = sqrt(r*N)*[0 0; 1 0];
    JN = chanconv(A,'kraus','choi');
end