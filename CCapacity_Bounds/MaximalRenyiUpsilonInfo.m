% compute the maximal Renyi Upsilon information
% input JN: the Choi matrix of the channel
%        d: dimension vector in the order of [da db]
%        l: divergence level l
% Writte by Kun Fang

function [t N0] = MaximalRenyiUpsilonInfo(JN,d,l)

da = d(1);
db = d(2);
JN = (JN+JN')/2; % avoid numerical problems from non-Hermitian JN

cvx_precision best

if l >= 1
cvx_begin sdp quiet

variable M(da*db,da*db) hermitian
variable N0(da*db,da*db) hermitian
variable N(da*db,da*db,l) hermitian
variable S(db,db) hermitian
variable R(da*db,da*db) hermitian
variable y

minimize y
subject to
    y*eye(da) - PartialTrace(M,2,[da db]) >= 0;
    S >= 0; R >= 0; y >= 0;
    [M JN; JN N(:,:,l)] >= 0;
    for i = 2:l
       [JN N(:,:,i); N(:,:,i) N(:,:,i-1)]  >= 0;
    end
    [JN N(:,:,1); N(:,:,1) N0]  >= 0;
    % set conditions
    R + PartialTranspose(N0,2,[da db]) >= 0;
    R - PartialTranspose(N0,2,[da db]) >= 0;
    Tensor(eye(da),S) + PartialTranspose(R,2,[da db]) >= 0;
    Tensor(eye(da),S) - PartialTranspose(R,2,[da db]) >= 0;    
    1 - trace(S) >= 0;
    
cvx_end
end

if l == 0
cvx_begin sdp quiet

variable M(da*db,da*db) hermitian
variable N0(da*db,da*db) hermitian
variable S(db,db) hermitian
variable R(da*db,da*db) hermitian
variable y

minimize y
subject to
    S >= 0; R >= 0; y >= 0;
    y*eye(da) - PartialTrace(M,2,[da db]) >= 0;
    [M JN; JN N0] >= 0;
    % set conditions
    R + PartialTranspose(N0,2,[da db]) >= 0;
    R - PartialTranspose(N0,2,[da db]) >= 0;
    Tensor(eye(da),S) + PartialTranspose(R,2,[da db]) >= 0;
    Tensor(eye(da),S) - PartialTranspose(R,2,[da db]) >= 0;
    1 - trace(S) >= 0;
cvx_end  
end
t = 2^l*log2(y);
end