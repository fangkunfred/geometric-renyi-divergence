% compute the maximal renyi Rains information 
% input JN: the Choi matrix of the channel
%        d: dimension vector in the order of [da db]
%        l: divergence level l
% Writte by Kun Fang


function t = MaximalRenyiRainsInfo(JN,d,l)

da = d(1);
db = d(2);
JN = (JN+JN')/2; % avoid numerical problems from non-Hermitian JN

cvx_precision best

if l >= 1
cvx_begin sdp quiet

variable K(da*db,da*db) complex
variable Z0(da*db,da*db) complex
variable Z(da*db,da*db,l) complex
variable W(da*db,da*db,l) hermitian
variable rho(da,da) hermitian

t = trace((K+K'-sum(W,3))*JN);

maximize t
subject to
    [Tensor(rho,eye(db)) K; K' Z(:,:,l)+Z(:,:,l)'] >= 0;
    for i = 2:l
        [W(:,:,i) Z(:,:,i); Z(:,:,i)' Z(:,:,i-1)+Z(:,:,i-1)'] >= 0;
    end
    [W(:,:,1) Z(:,:,1); Z(:,:,1)' Z0+Z0'] >= 0;
    Tensor(rho,eye(db)) + PartialTranspose(Z0+Z0',2,[da db]) >= 0;
    Tensor(rho,eye(db)) - PartialTranspose(Z0+Z0',2,[da db]) >= 0;
    trace(rho) - 1 == 0;

cvx_end
end

if l == 0
cvx_begin sdp quiet

variable K(da*db,da*db) complex
variable Z0(da*db,da*db) complex
variable rho(da,da) hermitian

t = trace((K+K')*JN);

maximize t
subject to
    [Tensor(rho,eye(db)) K; K' Z0+Z0'] >= 0;
    Tensor(rho,eye(db)) + PartialTranspose(Z0 + Z0',2,[da db]) >= 0;
    Tensor(rho,eye(db)) - PartialTranspose(Z0 + Z0',2,[da db]) >= 0;
    trace(rho) - 1 == 0;

cvx_end 
end
t = l*2^l-(2^l+1)*log2(2^l+1)+(2^l+1)*log2(t);
end



