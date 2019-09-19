% compute the max Upsilon information 
% input JN: the Choi matrix of the channel
%        d: dimension vector in the order of [da db]
% Writte by Kun Fang

function t = MaxUpsilonInfo(JN,d)

da = d(1);
db = d(2);
JN = (JN+JN')/2; % avoid numerical problems from non-Hermitian JN 

cvx_begin sdp quiet

variable S(db,db) hermitian
variable K(da*db,da*db) hermitian
variable R(da*db,da*db) hermitian

t = trace(S);

minimize t
subject to
    S >= 0; K >= 0; R >= 0;
    K - JN >= 0;
    R + PartialTranspose(K,2,[da db]) >= 0;
    R - PartialTranspose(K,2,[da db]) >= 0;
    Tensor(eye(da),S) + PartialTranspose(R,2,[da db]) >= 0;
    Tensor(eye(da),S) - PartialTranspose(R,2,[da db]) >= 0;
    
cvx_end
t = log2(t);
end