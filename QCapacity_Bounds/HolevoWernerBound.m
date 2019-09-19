% compute the holevo werner bound
% input JN: the Choi matrix of the channel
%        d: dimension vector in the order of [da db]
% Writte by Kun Fang

function r = HolevoWernerBound(JN,d)
da = d(1);
db = d(2);
JN = (JN+JN')/2; % avoid numerical problems from non-Hermitian JN

cvx_begin sdp quiet

variable X(da*db,da*db) complex
variable rho0(da,da) hermitian
variable rho1(da,da) hermitian

t = 1/2*trace(PartialTranspose(JN,2,[da db])*(X+X'));

maximize t
subject to
    [Tensor(rho0,eye(db)) X; X' Tensor(rho1,eye(db))] >= 0;
    rho0 >= 0; rho1 >= 0;
    trace(rho0) == 1; trace(rho1) == 1;
cvx_end
r = log2(t);
end