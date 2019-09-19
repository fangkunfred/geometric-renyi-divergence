% compute the max Rains information
% input JN: the Choi matrix of the channel
%        d: dimension vector in the order of [da db]
% Writte by Kun Fang

function r = MaxRainsInfo(JN,d)
da = d(1);
db = d(2);
JN = (JN+JN')/2; % avoid numerical problems from non-Hermitian JN

cvx_begin sdp quiet

variable R(da*db,da*db) hermitian
variable rho(da,da) hermitian

r = trace(JN*R);

maximize r
subject to
    R >= 0; rho >= 0; trace(rho) == 1;
    Tensor(rho,eye(db)) + PartialTranspose(R,2,[da db]) >= 0; 
    Tensor(rho,eye(db)) - PartialTranspose(R,2,[da db]) >= 0; 
           
cvx_end
r = log2(r);
end