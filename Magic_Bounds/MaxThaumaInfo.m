% compute the max Thauma information 
% input JN: the Choi matrix of the channel
%        d: dimension 
% Writte by Kun Fang

function t = MaxThaumaInfo(JN,d)

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
cvx_begin sdp quiet

variable V(da*db,da*db) hermitian
variable y

for a1 = 1:d
for a2 = 1:d
sum = 0;
for a3 = 1:d
for a4 = 1:d
sum = sum + abs(trace(Tensor(A(:,:,a1,a2),A(:,:,a3,a4))*V/d));
end
end
SUM(a1,a2)=(sum+sum')/2; % ensure that it is a real value
end 
end
         
minimize y
subject to
JN <= V;
% useless set 
for a1 = 1:d
for a2 = 1:d
SUM(a1,a2) <= y;
end
end
cvx_end

t = log2(y);
end