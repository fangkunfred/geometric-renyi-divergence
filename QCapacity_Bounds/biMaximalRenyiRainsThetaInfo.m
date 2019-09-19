% calculating the SDP relaxation of the bidirectional maximal renyi Rains Theta information 
% of a given quantum channel
% input quantum channel Choi matrix JN bipartite dimension [da db] hierarchy level l
% Writte by Kun Fang

function t = biMaximalRenyiRainsThetaInfo(JN,d,l)

dap = d(1);
dbp = d(2);
da = d(3);
db = d(4);
dtol = dap*dbp*da*db;
JN = (JN+JN')/2; % in case input JN is not strictly hermitian

if l >= 1
cvx_precision best
cvx_begin sdp quiet
   variable M(dtol,dtol) hermitian
   variable N0(dtol,dtol) hermitian
   variable N(dtol,dtol,l) hermitian
   variable R(dtol,dtol) hermitian
   variable y
   minimize y
   subject to
    PartialTrace(M,[3 4],[dap dbp da db]) <= y*eye(dap*dbp);
    [M JN; JN N(:,:,l)] >= 0;
    for i = 2:l
       [JN N(:,:,i); N(:,:,i) N(:,:,i-1)] >= 0;
    end
       [JN N(:,:,1); N(:,:,1) N0] >= 0;
       % theta set condition
       R + PartialTranspose(N0,[2 4],[dap dbp da db]) >= 0;
       R - PartialTranspose(N0,[2 4],[dap dbp da db]) >= 0;
       PartialTrace(R,[3 4],[dap dbp da db]) <= eye(dap*dbp);       
cvx_end
end

if l == 0
cvx_begin sdp quiet
variable M(dtol,dtol) hermitian
variable N0(dtol,dtol) hermitian
variable R(dtol,dtol) hermitian
variable y
minimize y
subject to
   PartialTrace(M,[3 4],[dap dbp da db]) <= y*eye(dap*dbp);
   [M JN; JN N0] >= 0;
   % theta set condition
   R + PartialTranspose(N0,[2 4],[dap dbp da db]) >= 0;
   R - PartialTranspose(N0,[2 4],[dap dbp da db]) >= 0;
   PartialTrace(R,[3 4],[dap dbp da db]) <= eye(dap*dbp);       
cvx_end
end
t = 2^l*log2(y);
end
