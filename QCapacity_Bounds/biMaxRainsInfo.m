% For given choi matrix of quantum channel N and its input output dimension
% compute its bidirectional max-Rains information bound
% the input dimension order is A'B'AB
% Writte by Kun Fang

function r = biMaxRainsInfo(JN,d)
dap = d(1);
dbp = d(2);
da = d(3);
db = d(4);
JN = (JN+JN')/2; % in case input JN is not strictly hermitian


cvx_begin sdp quiet
    variable V(dap*dbp*da*db,dap*dbp*da*db) hermitian
    variable Y(dap*dbp*da*db,dap*dbp*da*db) hermitian
    variable y
    minimize y
    subject to
        V >= 0; Y >= 0;
        PartialTranspose(V-Y,[2 4],[dap dbp da db]) >= JN;
        PartialTrace(V+Y,[3 4],[dap dbp da db]) <= y*eye(dap*dbp);
cvx_end
r = log2(y);
end