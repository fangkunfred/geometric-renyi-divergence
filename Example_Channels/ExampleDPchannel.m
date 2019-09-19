% Choi matrix of the depolarizing channel with parameter p and dimension d
function J = ExampleDPchannel(d,p)
J = (1-p)*MaxEntangled(d,0,0)*MaxEntangled(d,0,0)' + p*eye(d^2)/d;
end