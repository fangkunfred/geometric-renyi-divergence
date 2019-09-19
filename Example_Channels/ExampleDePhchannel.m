% Choi matrix of the dephasing channel with parameter p
function J = ExampleDePhchannel(p)
Phi = MaxEntangled(2,0,0)*MaxEntangled(2,0,0)';
    J = (1-p)*Phi+...
        p*Tensor(eye(2),Pauli('Z'))*Phi*Tensor(eye(2),Pauli('Z'));
end