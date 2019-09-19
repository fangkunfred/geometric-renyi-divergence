% compute the maximal renyi Thauma information 
% input JN: the Choi matrix of the channel
%        d: dimension 
%        l: divergence level l
% Writte by Kun Fang

function t = MaximalRenyiThaumaInfo(JN,d,l)

JN = (JN+JN')/2; % avoid numerical problems from non-Hermitian JN

% discrete wigner basis setup
X = GenPauli(1,0,d);
Z = GenPauli(0,1,d);
w = exp(2*pi*1j/d);
tau = exp((d+1)*pi*1j/d);

A0 = zeros(d);
for a1 = 1:d
for a2 = 1:d
T(:,:,a1,a2) = tau^(-a1*a2)*Z^(a1)*X^(a2);
A0 = A0 + T(:,:,a1,a2);
end
end
A0 = A0/d;

for a1 = 1:d
for a2 = 1:d
A(:,:,a1,a2) = T(:,:,a1,a2)*A0*T(:,:,a1,a2)';
end
end

da = d;
db = d;

cvx_precision best

if l >= 1
cvx_begin sdp quiet

variable M(da*db,da*db) hermitian
variable N0(da*db,da*db) hermitian
variable N(da*db,da*db,l) hermitian
variable y

for a1 = 1:d
for a2 = 1:d
sum = 0;
for a3 = 1:d
for a4 = 1:d
sum = sum + abs(trace(Tensor(A(:,:,a1,a2),A(:,:,a3,a4))*N0/d));
end
end
SUM(a1,a2)=(sum+sum')/2; % ensure that it is a real value
end 
end
         
minimize y
subject to
PartialTrace(M,2,[da db]) <= y*eye(da);
[M JN; JN N(:,:,l)] >= 0;
for i = 2:l
[JN N(:,:,i); N(:,:,i) N(:,:,i-1)] >= 0;
end
[JN N(:,:,1); N(:,:,1) N0] >= 0;
% useless set N0
for a1 = 1:d
for a2 = 1:d
SUM(a1,a2) <= 1;
end
end
cvx_end

end

if l == 0
cvx_begin sdp quiet

variable M(da*db,da*db) hermitian
variable N0(da*db,da*db) hermitian
variable y

for a1 = 1:d
for a2 = 1:d
sum = 0;
for a3 = 1:d
for a4 = 1:d
sum = sum + abs(trace(Tensor(A(:,:,a1,a2),A(:,:,a3,a4))*N0/d));
end
end
SUM(a1,a2)=(sum+sum')/2; % ensure that it is a real value
end 
end

minimize y
subject to
PartialTrace(M,2,[da db]) <= y*eye(da);
[M JN; JN N0] >= 0;
% useless set N0
for a1 = 1:d
for a2 = 1:d
SUM(a1,a2) <= 1;
end
end
cvx_end
end

t = 2^l*log2(y);
end