% Choi matrix of the amplitude damping channel with parameter gamma

function J = ExampleADchannel(gamma)
E{1} = [1 0; 0 sqrt(1-gamma)];
E{2} = [0 sqrt(gamma); 0 0];
J = chanconv(E,'kraus','choi');
end