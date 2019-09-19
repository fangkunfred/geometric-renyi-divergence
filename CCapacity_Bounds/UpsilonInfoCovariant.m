% compute the Upsilon information for covariant quantum channel
% input JN: the Choi matrix of the channel
%        d: dimension vector in the order of [da db]
% Writte by Kun Fang

function r = UpsilonInfoCovariant(JN,d)
da = d(1);
db = d(2);
JN = (JN+JN')/2; % avoid numerical problems from non-Hermitian JN

cvx_begin sdp quiet

variable JM(da*db,da*db) hermitian
variable R(da*db,da*db) hermitian
variable S(db,db) hermitian

r = quantum_rel_entr(JN/da,JM/da)/log(2);

minimize r
subject to
    JM >= 0; R >= 0; S >= 0;
    R + PartialTranspose(JM,2,[da db]) >= 0;
    R - PartialTranspose(JM,2,[da db]) >= 0;
    Tensor(eye(da),S) + PartialTranspose(R,2,[da db]) >= 0;
    Tensor(eye(da),S) - PartialTranspose(R,2,[da db]) >= 0;
    trace(S) <= 1;
    
cvx_end
end