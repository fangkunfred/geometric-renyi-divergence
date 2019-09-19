% For given choi matrix of quantum channel N and its input output dimension
% compute its bidirectional max-Rains Theta information bound
% Writte by Kun Fang

function r = biMaxRainsThetaInfo(JN,d)
dap = d(1);
dbp = d(2);
da = d(3);
db = d(4);
JN = (JN+JN')/2; % in case input JN is not strictly hermitian
cvx_begin sdp quiet
    variable R(dap*dbp*da*db,dap*dbp*da*db) hermitian
    variable K(dap*dbp*da*db,dap*dbp*da*db) hermitian   
    variable t    
    minimize t
    subject to
        JN <= K;
        R + PartialTranspose(K,[2 4],[dap dbp da db]) >= 0;
        R - PartialTranspose(K,[2 4],[dap dbp da db]) >= 0;
        PartialTrace(R,[3 4],[dap dbp da db]) <= t*eye(dap*dbp);
cvx_end
r = log2(t);
end